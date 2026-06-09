#!/usr/bin/env Rscript

#### Config ####
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: minitax.complete.R <config.tsv>", call. = FALSE)
}
config.file <- args[1]
#### ####
##

version <- 0.9

message("Welcome to", "\t", "𝑚ＩＮＩＴΛΧ", " !", "\n", "minitax, version: ", version)


message(
"
╔═════════════════════════════════╗
║                                 ║
║  ┌┬┐  ┬  ┐ ┌  ┬  ─┬─  ┌─┐  ╲╱   ║
║  │││  │  │┝│  │   │   │ │  ||   ║
║  ┴ ┴  ┴  ┘ └  ┴   ┴   ┴ ┴  ╱╲   ║
║                                 ║
╚═════════════════════════════════╝
")


#### Load packages quietely #####

suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(GenomicAlignments)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringi)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(future.apply)))

HAS_PHYLOSEQ <- requireNamespace("phyloseq", quietly = TRUE)
if (HAS_PHYLOSEQ) {
  suppressMessages(suppressWarnings(library(phyloseq)))
} else {
  message("phyloseq is not installed; phyloseq .rds outputs will be skipped.")
}

#### ####


#### ####
#config <- 'minitax_config_allNCBI.txt'   ##  -->> this is for debugging
config <- read.delim(config.file, stringsAsFactors = FALSE)
#print(config[,c(1:2)])

config <- config[,c(1:2)] %>% spread(argument, value)
t(config)

script.path <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
script.dir <- if (!is.na(script.path) && nzchar(script.path)) {
  dirname(normalizePath(script.path, mustWork = FALSE))
} else {
  getwd()
}

is_config_na <- function(value) {
  length(value) == 0 || is.na(value) || !nzchar(trimws(value)) || toupper(trimws(value)) == "NA"
}

set_config_default <- function(name, default) {
  if (!name %in% names(config) || is_config_na(config[[name]][1])) {
    config[[name]] <<- default
  }
}

split_config_values <- function(value, default = character()) {
  if (is_config_na(value)) return(default)
  values <- trimws(unlist(strsplit(value, "[;,]")))
  values[nzchar(values)]
}

parse_logical_config <- function(value, default = FALSE) {
  if (is_config_na(value)) return(default)
  value <- toupper(trimws(value))
  if (value %in% c("TRUE", "T", "YES", "Y", "1")) return(TRUE)
  if (value %in% c("FALSE", "F", "NO", "N", "0")) return(FALSE)
  stop("Invalid logical config value: ", value, call. = FALSE)
}

required.config <- c("platform", "db", "db.dir", "nproc")
missing.config <- required.config[!required.config %in% names(config)]
present.required <- intersect(required.config, names(config))
missing.config <- unique(c(
  missing.config,
  present.required[vapply(config[present.required], is_config_na, logical(1))]
))
if (length(missing.config) > 0) {
  stop("Missing required config argument(s): ", paste(missing.config, collapse = ", "), call. = FALSE)
}

set_config_default("minitax.dir", script.dir)
set_config_default("misc.dir", script.dir)
set_config_default("project", "project")
set_config_default("Vregion", "WGS")
set_config_default("outdir", "NA")
set_config_default("outputs", "bam.sum")
set_config_default("methods", "BestAln")
set_config_default("mapq.filt", NA)
set_config_default("best.mapq", "TRUE")
set_config_default("keep.max.cigar", "TRUE")
set_config_default("keep.highest.mapq.aln.only", "TRUE")
set_config_default("crop.na.tax", "FALSE")
set_config_default("multicore", "TRUE")
set_config_default("saveRAM", "FALSE")
set_config_default("reuse.taxa.cache", "TRUE")
set_config_default("dt.threads", "1")
set_config_default("CIGAR_points", "match_score = 1; mismatch_score = -3; insertion_score = -2; deletion_score = -2; gap_opening_penalty = -2; gap_extension_penalty = -1")

nproc <- as.integer(config$nproc) #args[2] #32 #
if (is.na(nproc) || nproc < 1) {
  stop("Config argument nproc must be a positive integer.", call. = FALSE)
}
#### ####
#stop()

#### Convert paths, according to OS #####
convert_path <- function(path) {
  gsub("^/mnt/([a-zA-Z])/","\\U\\1:/", path, perl = TRUE)
}

if (Sys.info()[["sysname"]] == 'Windows') {
  for (path.arg in intersect(c("minitax.dir", "misc.dir", "db.dir"), names(config))) {
    config[[path.arg]] <- convert_path(config[[path.arg]])
  }
}
#### ####
##

#### functions ####
## minitax
minitax.dir <- config$minitax.dir ## 'my.R.packages/minitax'
source_required <- function(base.dir, rel.path) {
  path <- file.path(base.dir, rel.path)
  if (!file.exists(path)) {
    stop("Required source file not found: ", path, call. = FALSE)
  }
  source(path)
}

source_required(minitax.dir, 'R/minitax.v2.R')
source_required(minitax.dir, 'R/cigar_score.R')
source_required(minitax.dir, 'R/cigar_to_length.R')
#source(paste0( minitax.dir, '/R/PS_from_taxa.sum.R'))
#source(paste0( minitax.dir, '/R/minitax.R'))
#source(paste0( minitax.dir, '/R/tax.identity.ranks.R'))
#source(paste0( minitax.dir, '/R/make.krona.out.R'))
#source(paste0( minitax.dir, '/R/estimate.species.R'))
#source(paste0( minitax.dir, '/R/import.bam.R'))
source_required(minitax.dir, 'R/cigar.sum.R')
#source(paste0( minitax.dir, '/R/get_chunks.R'))
#source(paste0( minitax.dir, '/R/taxprofile.from.minitax.R'))
source_required(minitax.dir, 'R/rename_duptaxa.R')
source_required(minitax.dir, 'R/minitax.wrapfun.complete.R')
source_required(minitax.dir, 'R/BestAln.R')
source_required(minitax.dir, 'R/add_lineage.R')

## misc package
misc.dir <- config$misc.dir ## 'my.R.packages/Rlyeh-main'
source_required(misc.dir, 'R/ov.from.bam2.R')
#source(paste0( misc.dir, '/R/ov.from.bam3.R'))
source_required(misc.dir, 'R/dt.from.bam.R')
source_required(misc.dir, 'R/get.best.aln.R')
source_required(misc.dir, 'R/misc.functions.R')
source_required(misc.dir, 'R/add_rightmost_non.NA_col.R')
source_required(misc.dir, 'R/paste_columns_dt.R')
bam.flags <- read.delim(file.path(misc.dir, 'R/bam.flags.tsv'))
#### ####
##


#### Project options ####
db        <- config$db ##'proGcontigs'
db.dir    <- config$db.dir ## paste0(path.prefix, 'data/databases/proGenomes')
project   <- config$project ## 'MCM'
Vregion   <- config$Vregion ## 'v1_2'
platform  <- config$platform ##'illumina'
default.outdir <- paste('minitax', project, platform, Vregion, db, sep = '_')
outdir    <- if (is_config_na(config$outdir)) default.outdir else config$outdir
#outdir <- paste0(outdir, '.partII')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir.label <- basename(gsub("[/\\\\]+$", "", outdir))
if (is_config_na(outdir.label)) outdir.label <- "minitax"
method_result_path <- function(methods, suffix) {
  file.path(outdir, paste0(outdir.label, "_", methods, suffix))
}


bam.all.out      <- paste0(outdir, '/', 'bam.all.tsv')
bamstats.all.out <- paste0(outdir, '/', 'bamstats.tsv')
taxa.all.out     <- paste0(outdir, '/', 'minitax.tsv')
taxa.prog.out    <- paste0(outdir, '/', 'taxa.tsv')
taxa.comp.out    <- paste0(outdir, '/', 'all.reads.tsv')
tax.abund.out    <- paste0(outdir, '/', 'tax.abund.tsv')
spec.est.all.out <- paste0(outdir, '/', 'species.estimates.tsv')
#### ####
##

#### Database import ####
if (db == 'proGcontigs_2') {
  prog.db <- fread(paste0(db.dir, '/proGenomes2.1_specI_lineageNCBI.tab'), header = F)
  colnames(prog.db) <- c("genome", "superkingdom", "phylum", "class", "order", "family", "genus", "species")
  prog.db <- data.frame(taxid=gsub('\\..*', '', prog.db$genome  ), prog.db)
  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  prog.db$taxid <- as.integer(prog.db$taxid)

  prog.db[,ranks]  <- as.data.frame(apply(prog.db[,ranks], 2, function(x) str_replace(as.character(unlist(x)),  pattern = "^\\d+\\s", replacement = "")))

  prog.gt <- gather(prog.db, rank, taxon, -c(1:2))
  prog.db.spec.uni <- data.table(unique.data.frame(prog.db[,ranks]))

  db.data <- data.table(prog.db)
  db.name <- db #'proGcontigs_2'

  prog.db.uni     <- db.data %>% distinct(across(all_of(c('taxid', ranks)))) #unique.data.frame(db.data[,c('taxid', ranks)])
  db.uni.data     <- data.table(prog.db.uni)

} else if (db == 'EMUdb') {

  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  emu.db    <- read.delim(paste0(db.dir, '/taxonomy.tsv'))
  colnames(emu.db)[1] <- 'taxid'
  emu.db <- emu.db[,c('taxid', ranks)]
  emu.db[emu.db == ''] <- NA

  emu.fasta <- seqinr::read.fasta(paste0(db.dir, '/species_taxid.fasta'), whole.header = T)
  emu.fasta <- data.frame(names=names(emu.fasta), taxid=gsub(':.*', '', names(emu.fasta)))
  emu.fasta$seq_id <- gsub(' \\[.*', '', emu.fasta$names)

  emu.idx    <- merge(emu.fasta, emu.db, by='taxid')

  db.data <- data.table(emu.db)
  db.name <- db #'EMUdb'
  db.uni.data <- db.data

} else if (db == 'rrn') {
  NULL
} else if (db=='proGcontigs_3.host' | db=='proGcontigs_3.repres') {

  prog.db <- fread(paste0(db.dir, '/progenomes3.db.tsv'), header = T)
  ranks   <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  colnames(prog.db)[7:13] <- ranks
  #c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  #prog.db <- data.frame(taxid=gsub('\\..*', '', prog.db$genome  ), prog.db)
  #prog.db$taxid <- as.integer(prog.db$taxid)
  #prog.gt <- gather(prog.db, rank, taxon, -c(1:2))
  prog.db.spec.uni <- unique(prog.db[,..ranks])

  db.data <- prog.db
  db.name <- db

  prog.db.uni     <- db.data %>% distinct(across(all_of(c('seqnames', ranks)))) #unique.data.frame(db.data[,c('taxid', ranks)])
  db.uni.data     <- prog.db.uni

} else if (db=="all_NCBI_genomes") {

  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

  db.data     <- fread(paste0(db.dir, "/NCBI.db.tsv"), header = T)
  db.uni.data <- fread(paste0(db.dir, "/NCBI.db.uni.tsv"), header = T)
  ## we need taxid from db data!
  db.uni.data <- db.data
  db.name <- db

} else if (db == "GTDB_SSU") {
  
  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  
  gtdb.tax.path <- file.path(db.dir, "gtdb_ssu_taxonomy.tsv")
  
  if (!file.exists(gtdb.tax.path)) {
    stop("GTDB_SSU database selected, but taxonomy file was not found: ",
         gtdb.tax.path, call. = FALSE)
  }
  
  gtdb.tax <- fread(
    gtdb.tax.path,
    header = FALSE,
    sep = "\t",
    col.names = c("seqnames", "lineage_raw"),
    na.strings = c("", "NA")
  )
  
  # Split GTDB lineage:
  # d__Bacteria;p__Bacteroidota;c__Bacteroidia;...
  lineage_split <- tstrsplit(gtdb.tax$lineage_raw, ";", fixed = TRUE)
  
  if (length(lineage_split) < length(ranks)) {
    stop("GTDB lineage has fewer than ", length(ranks), " ranks in ",
         gtdb.tax.path, call. = FALSE)
  }
  
  for (i in seq_along(ranks)) {
    gtdb.tax[, (ranks[i]) := lineage_split[[i]]]
  }
  
  # Optional but recommended: remove GTDB rank prefixes d__, p__, c__, ...
  for (r in ranks) {
    gtdb.tax[, (r) := sub("^[a-z]__", "", get(r))]
    gtdb.tax[get(r) == "", (r) := NA_character_]
  }
  
  # taxid is not numeric here; use the GTDB genome accession as stable taxid-like ID
  gtdb.tax[, taxid := seqnames]
  
  db.data <- gtdb.tax[, c("seqnames", "taxid", ranks), with = FALSE]
  db.uni.data <- unique(db.data)
  db.name <- db
  
  setkeyv(db.uni.data, "seqnames")

} else {

  db.data     <- fread(paste0(db.dir, "/db_data.tsv"), header = T, na.strings = '')
  db.uni.data <- db.data
  db.name     <- db

  message('If another database was used, then there should be a "db_data.tsv" file present in the
          database directory, containing the taxonomies (lineages) of the sequences of the databases (the header of the fasta files).
          In this file, the two column names (first row) should be:',
          'seqnames', 'taxid', '\n',
          'That is, the first coulmn is the sequence identifers of the dabase fastafile,
           the second is the identifier of the organism,
           and the followings are the taxonomies.')

  try({ ranks <- unlist(strsplit(config$ranks, ', ')) })
  if (!exists('ranks')) { ranks <- colnames(db.uni.data)[-c(1)] }
  # Identify the last entry of the ranks vector
  last_rank <- tail(ranks, n = 1)

  #c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  cols  <- c('seqnames', ranks) ## 'taxid',


}


message(db.name, ' database loaded')

message('The following taxon levels will be used: ', paste0(ranks, collapse='; '))

if (is.data.table(db.uni.data)) {
  if ('seqnames' %in% colnames(db.uni.data)) {
    setkeyv(db.uni.data, 'seqnames')
  } else if ('taxid' %in% colnames(db.uni.data)) {
    setkeyv(db.uni.data, 'taxid')
  }
}



#### ####
##


#### Metadata and parent directory of bamfiles from minimap2 output ####
pardir    <- file.path(outdir, 'bam')
pattern   <- '.bam'
bamfiles  <- if (dir.exists(pardir)) {
  grep('\\.bai$', list.files(pardir, pattern = pattern), value = T, invert = T)
} else {
  character()
}
cached.taxa.files <- list.files(
  file.path(outdir, 'best_alignments_w_taxa'),
  '.*_best_alignments_w_taxa\\.tsv$',
  full.names = TRUE
)

if (length(bamfiles) == 0 && length(cached.taxa.files) == 0) {
  stop("No BAM files found in ", pardir, " and no cached best_alignments_w_taxa files found.", call. = FALSE)
}

### create results directories
outputs <- split_config_values(config$outputs, default = 'bam.sum')
for(resdir in outputs) {
  dir.create(file.path(outdir, resdir), recursive = TRUE, showWarnings = FALSE)
}

## filter
#bamfiles <- bamfiles[14]

if (length(bamfiles) > 0) {
  samples <- sub('\\.bam$', '', basename(bamfiles))
  metadata.bamfiles <- bamfiles
} else {
  samples <- gsub('_best_alignments_w_taxa.tsv$', '', basename(cached.taxa.files))
  metadata.bamfiles <- NA
}

#source('make.metadata.WGS.R')
#metadata$db <- db

metadata <- data.frame(sample=samples, bamfile=metadata.bamfiles,
                       project=project,
                       db=db,
                       Vregion=Vregion,
                       platform=platform,
                       row.names = samples,
                       workflow = 'minitax',
                       mapq.filt =  config$mapq.filt
)

write_tsv(metadata, paste0(outdir, "/metadata.tsv"))
#metadata <- read.delim(paste0(outdir, "/metadata.tsv"))
#### ####
##


#### Minitax settings
keep.highest.mapq.aln.only <- parse_logical_config(config$keep.highest.mapq.aln.only, default = TRUE)
crop.na.tax <- parse_logical_config(config$crop.na.tax, default = FALSE)
mapq.filt   <- config$mapq.filt; if(!is_config_na(mapq.filt)) {
  mapq.start <- as.integer(gsub(':.*', '', config$mapq.filt))
  mapq.end <- as.integer(gsub('.*:', '', config$mapq.filt))
  if (is.na(mapq.start) || is.na(mapq.end)) {
    stop("mapq.filt must be NA, a single integer, or an integer range like 0:59.", call. = FALSE)
  }
  mapq.filt   <- mapq.start:mapq.end  # default: NA
} else {
  mapq.filt <- NA
}
methods.to.use <- split_config_values(config$methods, default = 'BestAln')
best.mapq      <- parse_logical_config(config$best.mapq, default = TRUE)
multicore      <- parse_logical_config(config$multicore, default = TRUE)
saveRAM        <- parse_logical_config(config$saveRAM, default = FALSE)
reuse.taxa.cache <- parse_logical_config(config$reuse.taxa.cache, default = TRUE)
dt.threads     <- as.integer(config$dt.threads)
if (is.na(dt.threads) || dt.threads < 1) {
  stop("Config argument dt.threads must be a positive integer.", call. = FALSE)
}
data.table::setDTthreads(dt.threads)

keep.max.cigar <- parse_logical_config(config$keep.max.cigar, default = TRUE)
CIGAR_points   <- separate(data.frame(stringr::str_split_1(config$CIGAR_points, ';\\s*')), col = 1, into = c('key', 'value'), sep = '\\s*=\\s*')
CIGAR_points$value <- as.integer(CIGAR_points$value)
if (any(is.na(CIGAR_points$value))) {
  stop("CIGAR_points contains non-integer values.", call. = FALSE)
}
CIGAR_points <- CIGAR_points %>% spread(key, value)

runtime.outputs <- outputs
if (length(bamfiles) > 0 || length(methods.to.use) > 1) {
  runtime.outputs <- unique(c(runtime.outputs, 'best_alignments_w_taxa'))
}
for (resdir in runtime.outputs) {
  dir.create(file.path(outdir, resdir), recursive = TRUE, showWarnings = FALSE)
}


#### configs used
message('Configs used: ')
print(t(config))
##

####
message('Constructed metadata:')
print(metadata)

#### Set up future to use multiple cores with future
configure_minitax_plan <- function(task_count) {
  workers <- max(1L, min(nproc, max(1L, task_count)))
  if (!multicore) {
    future::plan(future::sequential)
    message("future plan: sequential")
    return(1L)
  }
  if (Sys.info()[["sysname"]] == 'Windows') {
    future::plan(future::multisession, workers = workers)
    message("future plan: multisession with ", workers, " worker(s)")
  } else {
    future::plan(future::multicore, workers = workers)
    message("future plan: multicore with ", workers, " worker(s)")
  }
  workers
}

timing.records <- list()
record_timing <- function(timing.dt) {
  if (is.null(timing.dt) || nrow(timing.dt) == 0) return(invisible(NULL))
  timing.records[[length(timing.records) + 1L]] <<- timing.dt
  fwrite(rbindlist(timing.records, fill = TRUE),
         file.path(outdir, 'classification_timing.tsv'),
         sep = '\t', na = 'NA')
  invisible(NULL)
}

write_method_outputs <- function(methods, taxa.sum) {
  setDT(taxa.sum)
  saveRDS(taxa.sum, method_result_path(methods, '_taxa.all.sum.rds'))
  fwrite(taxa.sum, method_result_path(methods, '_taxa.all.sum.tsv'), sep = '\t', na = 'NA')

  message('output of ', methods, ': ')
  print(head(taxa.sum))

  if (!HAS_PHYLOSEQ) {
    message("Skipping phyloseq object for ", methods, " because phyloseq is not installed.")
    return(invisible(NULL))
  }

  ps <- NULL
  taxa.for.ps <- copy(taxa.sum)
  try({
    if (methods == 'SpeciesEstimate') {
      taxa.for.ps[, count := norm_count]
    }

    sample_data        <- metadata
    sample_data$method <- methods

    taxa.for.ps <- taxa.for.ps[, c('tax.identity', 'lineage', ranks, 'sample', 'count'), with = FALSE]
    taxa.for.ps[is.na(tax.identity), c('tax.identity', 'lineage', ranks) := 'unclassified']

    otutab <- dcast.data.table(taxa.for.ps, lineage ~ sample, value.var = 'count', fill = 0)
    otutab <- data.frame(otutab[, -1, with = FALSE], row.names = otutab$lineage)

    taxtab <- data.frame(unique(taxa.for.ps[, c('lineage', ranks, 'tax.identity'), with = FALSE]))
    rownames(taxtab) <- taxtab$lineage
    taxtab <- taxtab[rownames(otutab), ]

    ps <- phyloseq(otu_table(as.matrix(otutab), taxa_are_rows = TRUE),
                   tax_table(as.matrix(taxtab)),
                   sample_data(sample_data))

    sample_names(ps) <- paste(sample_data(ps)$workflow, sample_data(ps)$db,
                              sample_data(ps)$method, sample_names(ps), sep = '_')

    message('Generated phyloseq object: ')
    print(ps)
    saveRDS(ps, method_result_path(methods, '_PS.rds'))
  })
  if (is.null(ps)) stop("phyloseq object could not be created!")
  invisible(NULL)
}

run_minitax_method <- function(methods, taxa.files) {
  message('\n \n ', methods, ' method started, on: ', date())
  message('Running ', methods, ' from cached tax-annotated alignments: \n',
          paste(taxa.files, collapse = '\n'))

  configure_minitax_plan(length(taxa.files))
  method.results <- future_lapply(
    taxa.files,
    summarise_minitax_taxa_file,
    config = config,
    methods = methods,
    keep.max.cigar = keep.max.cigar,
    CIGAR_points = CIGAR_points,
    mapq.filt = mapq.filt,
    outputs = outputs,
    best.mapq = best.mapq,
    db = db,
    db.data = db.data,
    db.uni.data = db.uni.data,
    ranks = ranks,
    outdir = outdir,
    dt.threads = dt.threads,
    future.seed = TRUE
  )

  record_timing(rbindlist(lapply(method.results, `[[`, 'timing'), fill = TRUE))
  taxa.sum <- rbindlist(lapply(method.results, `[[`, 'taxa.sum'), fill = TRUE)

  if (methods == 'LCA') {
    taxa.sum <- rename_duptaxa(taxa.sum, ranks = ranks)
  } else if (methods == 'SpeciesEstimate') {
    taxID <- ranks[length(ranks)]
    taxa.sum[, tax.identity := get(taxID)]
    taxa.sum[is.na(tax.identity), tax.identity := 'unclassified']
    taxa.sum[, lineage := tax.identity]
  }

  write_method_outputs(methods, taxa.sum)
  if (saveRAM) {
    rm(method.results, taxa.sum)
    gc()
  }
  message(methods, ' method finished, on: ', date(), '\n \n ')
  invisible(NULL)
}

####
message('everything is ready! Starting minitax!')

#### RUN MINITAX ON ALIGNMENTS ####
taxa.files <- character()
if (length(bamfiles) > 0) {
  bam.paths <- file.path(pardir, bamfiles)
  message('Preparing reusable best_alignments_w_taxa cache from BAM files: \n',
          paste(bam.paths, collapse = '\n'))

  configure_minitax_plan(length(bam.paths))
  cache.results <- future_lapply(
    bam.paths,
    prepare_minitax_taxa_cache,
    config = config,
    keep.max.cigar = keep.max.cigar,
    CIGAR_points = CIGAR_points,
    mapq.filt = mapq.filt,
    outputs = runtime.outputs,
    best.mapq = best.mapq,
    db = db,
    db.data = db.data,
    db.uni.data = db.uni.data,
    ranks = ranks,
    outdir = outdir,
    reuse.taxa.cache = reuse.taxa.cache,
    dt.threads = dt.threads,
    future.seed = TRUE
  )

  record_timing(rbindlist(lapply(cache.results, `[[`, 'timing'), fill = TRUE))
  taxa.files <- vapply(cache.results, `[[`, character(1), 'cache_file')
  if (saveRAM) {
    rm(cache.results)
    gc()
  }
} else {
  taxa.files <- cached.taxa.files
  message('No BAM files found; using existing cached tax-annotated alignments.')
}

taxa.files <- taxa.files[file.exists(taxa.files)]
if (length(taxa.files) == 0) {
  stop("No best_alignments_w_taxa cache files are available for summarisation.", call. = FALSE)
}

for (methods in methods.to.use) {
  run_minitax_method(methods, taxa.files)
}

if (length(timing.records) > 0) {
  fwrite(rbindlist(timing.records, fill = TRUE),
         file.path(outdir, 'classification_timing.tsv'),
         sep = '\t', na = 'NA')
  message('Classification timing written to: ', file.path(outdir, 'classification_timing.tsv'))
}

future::plan(future::sequential)
#### ####

##############################
