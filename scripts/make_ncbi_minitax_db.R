#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

usage <- function() {
  cat(
    "Usage: make_ncbi_minitax_db.R --assemblies FILE --seq-map FILE --seq-lengths FILE --taxdump-dir DIR --outdir DIR\n",
    file = stderr()
  )
}

args <- commandArgs(trailingOnly = TRUE)
opts <- list()
i <- 1
while (i <= length(args)) {
  key <- args[[i]]
  if (!startsWith(key, "--") || i == length(args)) {
    usage()
    stop("Invalid arguments.", call. = FALSE)
  }
  opts[[sub("^--", "", key)]] <- args[[i + 1]]
  i <- i + 2
}

required <- c("assemblies", "seq-map", "seq-lengths", "taxdump-dir", "outdir")
missing <- required[!required %in% names(opts)]
if (length(missing) > 0) {
  usage()
  stop("Missing required option(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

assemblies_file <- opts[["assemblies"]]
seq_map_file <- opts[["seq-map"]]
seq_lengths_file <- opts[["seq-lengths"]]
taxdump_dir <- opts[["taxdump-dir"]]
outdir <- opts[["outdir"]]

for (path in c(assemblies_file, seq_map_file, seq_lengths_file, file.path(taxdump_dir, "nodes.dmp"), file.path(taxdump_dir, "names.dmp"))) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

message("Reading assembly and sequence metadata...")
assemblies <- fread(assemblies_file, sep = "\t", quote = "", na.strings = c("", "na", "NA"))
seq_map <- fread(seq_map_file, sep = "\t", quote = "", na.strings = c("", "na", "NA"))
seq_lengths <- fread(seq_lengths_file, sep = "\t", header = FALSE, col.names = c("seqnames", "seqlengths"))

setnames(seq_map, c("seqnames", "ident", "taxid"))
seq_map[, taxid := as.character(taxid)]
assemblies[, taxid := as.character(taxid)]
seq_lengths[, seqlengths := as.integer(seqlengths)]

message("Reading NCBI taxdump...")
nodes <- fread(
  file.path(taxdump_dir, "nodes.dmp"),
  sep = "|",
  header = FALSE,
  quote = "",
  strip.white = TRUE,
  select = 1:3,
  col.names = c("taxid", "parent_taxid", "rank")
)
names <- fread(
  file.path(taxdump_dir, "names.dmp"),
  sep = "|",
  header = FALSE,
  quote = "",
  strip.white = TRUE,
  select = c(1, 2, 4),
  col.names = c("taxid", "name_txt", "name_class")
)

nodes[, taxid := trimws(as.character(taxid))]
nodes[, parent_taxid := trimws(as.character(parent_taxid))]
nodes[, rank := trimws(as.character(rank))]
names[, taxid := trimws(as.character(taxid))]
names[, name_txt := trimws(as.character(name_txt))]
names[, name_class := trimws(as.character(name_class))]
names <- names[name_class == "scientific name", .(taxid, name_txt)]

parent_by_taxid <- setNames(nodes$parent_taxid, nodes$taxid)
rank_by_taxid <- setNames(nodes$rank, nodes$taxid)
name_by_taxid <- setNames(names$name_txt, names$taxid)

lineage_for_taxid <- function(taxid) {
  lineage <- setNames(rep(NA_character_, length(ranks)), ranks)
  current <- as.character(taxid)
  seen <- character()

  while (!is.na(current) && nzchar(current) && !(current %in% seen)) {
    seen <- c(seen, current)
    rank <- rank_by_taxid[[current]]
    if (!is.null(rank) && rank %in% ranks && is.na(lineage[[rank]])) {
      tax_name <- name_by_taxid[[current]]
      if (is.null(tax_name)) tax_name <- NA_character_
      lineage[[rank]] <- tax_name
    }
    parent <- parent_by_taxid[[current]]
    if (is.null(parent) || is.na(parent) || parent == current) break
    current <- parent
  }

  result <- as.data.table(as.list(lineage))
  result[, taxid := as.character(taxid)]
  setcolorder(result, c("taxid", ranks))
  result
}

message("Computing lineages...")
lineages <- rbindlist(lapply(unique(seq_map$taxid), lineage_for_taxid), fill = TRUE)

db <- merge(seq_map, lineages, by = "taxid", all.x = TRUE)
db <- merge(
  db,
  assemblies[, .(ident = assembly_accession, organism_name, tax_group)],
  by = "ident",
  all.x = TRUE
)
db[tax_group == "viral" & is.na(superkingdom), superkingdom := "Viruses"]
db[tax_group == "viral" & is.na(species) & !is.na(organism_name), species := organism_name]
db[, c("organism_name", "tax_group") := NULL]
setcolorder(db, c("seqnames", "ident", "taxid", ranks))
setkey(db, seqnames)

message("Writing NCBI.db.tsv...")
fwrite(db, file.path(outdir, "NCBI.db.tsv"), sep = "\t", quote = FALSE, na = "")

db_uni <- unique(db[, c("seqnames", "taxid", ranks), with = FALSE])
fwrite(db_uni, file.path(outdir, "NCBI.db.uni.tsv"), sep = "\t", quote = FALSE, na = "")

db_uni_spec <- unique(db[, ranks, with = FALSE])
setorder(db_uni_spec, superkingdom, phylum, class, order, family, genus, species)
fwrite(db_uni_spec, file.path(outdir, "NCBI.db.uni.spec.tsv"), sep = "\t", quote = FALSE, na = "")

message("Computing genome sizes...")
seq_len_map <- merge(seq_map, seq_lengths, by = "seqnames", all.x = TRUE)
genome_sizes <- seq_len_map[, .(genome_size = sum(seqlengths, na.rm = TRUE)), by = .(ident, taxid)]

assembly_info <- unique(assemblies[, .(
  ident = assembly_accession,
  organism_name,
  infraspecific_name,
  isolate,
  assembly_level,
  source,
  tax_group
)])
assembly_info[, infraspecific_name := as.character(infraspecific_name)]
assembly_info[, isolate := as.character(isolate)]
assembly_info[, strain := fifelse(
  !is.na(infraspecific_name) & nzchar(infraspecific_name),
  sub("^strain=", "", infraspecific_name),
  fifelse(!is.na(isolate) & nzchar(isolate), isolate, "")
)]

genome_sizes <- merge(
  genome_sizes,
  assembly_info[, .(ident, strain, organism_name, tax_group)],
  by = "ident",
  all.x = TRUE
)
genome_sizes <- merge(genome_sizes, lineages, by = "taxid", all.x = TRUE)
genome_sizes[tax_group == "viral" & is.na(superkingdom), superkingdom := "Viruses"]
genome_sizes[tax_group == "viral" & is.na(species) & !is.na(organism_name), species := organism_name]
genome_sizes[, c("organism_name", "tax_group") := NULL]
setcolorder(genome_sizes, c(ranks, "strain", "taxid", "ident", "genome_size"))
setorder(genome_sizes, superkingdom, phylum, class, order, family, genus, species, ident)
fwrite(genome_sizes, file.path(outdir, "NCBI.db.genomesize.tsv"), sep = "\t", quote = FALSE, na = "")

message("Wrote:")
message("  ", file.path(outdir, "NCBI.db.tsv"))
message("  ", file.path(outdir, "NCBI.db.uni.tsv"))
message("  ", file.path(outdir, "NCBI.db.uni.spec.tsv"))
message("  ", file.path(outdir, "NCBI.db.genomesize.tsv"))
