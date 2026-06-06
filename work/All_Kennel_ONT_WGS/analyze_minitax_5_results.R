suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% ""
script_path <- sub("^--file=", "", file_arg)
project_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
if (!dir.exists(project_dir) || basename(project_dir) != "All_Kennel_ONT_WGS") {
  project_dir <- getwd()
}
setwd(project_dir)

outdir <- "minitax_bacteria_refseq_5"
analysis_dir <- file.path(outdir, "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

method_files <- list(
  BestAln = file.path(outdir, paste0(outdir, "_BestAln_taxa.all.sum.rds")),
  SpeciesEstimate = file.path(outdir, paste0(outdir, "_SpeciesEstimate_taxa.all.sum.rds"))
)

read_method <- function(method_name, path) {
  if (!file.exists(path)) {
    stop("Missing result file for ", method_name, ": ", path, call. = FALSE)
  }
  dt <- as.data.table(readRDS(path))
  dt[, method := method_name]
  if (!"count" %in% names(dt)) dt[, count := NA_real_]
  if (!"norm_count" %in% names(dt)) dt[, norm_count := NA_real_]
  dt[, abundance := fifelse(!is.na(count), as.numeric(count), as.numeric(norm_count))]
  dt[, species_label := fifelse(!is.na(species) & species != "", species,
                                fifelse(!is.na(tax.identity) & tax.identity != "", tax.identity, "unclassified"))]
  dt[, genus_label := fifelse(!is.na(genus) & genus != "", genus, "unclassified")]
  dt[, phylum_label := fifelse(!is.na(phylum) & phylum != "", phylum, "unclassified")]
  dt[]
}

results <- rbindlist(Map(read_method, names(method_files), method_files), fill = TRUE)

for (method_name in names(method_files)) {
  dt <- results[method == method_name]
  export_dt <- copy(dt)
  export_dt[, method := NULL]
  fwrite(export_dt, file.path(outdir, paste0(outdir, "_", method_name, "_taxa.all.sum.tsv")), sep = "\t", na = "NA")
}

sample_totals <- results[, .(
  total_abundance = sum(abundance, na.rm = TRUE),
  taxa_rows = .N,
  species_detected = uniqueN(species_label[species_label != "unclassified"])
), by = .(method, sample)]
sample_totals[, total_abundance := round(total_abundance, 3)]
fwrite(sample_totals, file.path(analysis_dir, "sample_totals_by_method.tsv"), sep = "\t", na = "NA")

species_sum <- results[, .(abundance = sum(abundance, na.rm = TRUE)),
                       by = .(method, sample, phylum_label, genus_label, species_label)]
species_sum[, total := sum(abundance), by = .(method, sample)]
species_sum[, rel_abundance_pct := fifelse(total > 0, 100 * abundance / total, NA_real_)]
setorder(species_sum, method, sample, -rel_abundance_pct, species_label)
species_sum[, rank := seq_len(.N), by = .(method, sample)]
top_species <- species_sum[rank <= 20]
top_species[, `:=`(
  abundance = round(abundance, 3),
  total = round(total, 3),
  rel_abundance_pct = round(rel_abundance_pct, 3)
)]
fwrite(top_species, file.path(analysis_dir, "top20_species_by_method.tsv"), sep = "\t", na = "NA")

genus_sum <- results[, .(abundance = sum(abundance, na.rm = TRUE)),
                     by = .(method, sample, phylum_label, genus_label)]
genus_sum[, total := sum(abundance), by = .(method, sample)]
genus_sum[, rel_abundance_pct := fifelse(total > 0, 100 * abundance / total, NA_real_)]
setorder(genus_sum, method, sample, -rel_abundance_pct, genus_label)
genus_sum[, rank := seq_len(.N), by = .(method, sample)]
top_genera <- genus_sum[rank <= 20]
top_genera[, `:=`(
  abundance = round(abundance, 3),
  total = round(total, 3),
  rel_abundance_pct = round(rel_abundance_pct, 3)
)]
fwrite(top_genera, file.path(analysis_dir, "top20_genera_by_method.tsv"), sep = "\t", na = "NA")

comparison <- dcast(
  species_sum,
  sample + species_label ~ method,
  value.var = "rel_abundance_pct",
  fill = 0
)
if (all(c("BestAln", "SpeciesEstimate") %in% names(comparison))) {
  comparison[, delta_speciesestimate_minus_bestaln := SpeciesEstimate - BestAln]
  comparison[, abs_delta := abs(delta_speciesestimate_minus_bestaln)]
  setorder(comparison, sample, -abs_delta)
}
fwrite(comparison, file.path(analysis_dir, "species_method_relative_abundance_comparison.tsv"), sep = "\t", na = "NA")

fmt_top <- function(method_name, sample_name, n = 8) {
  dt <- top_species[method == method_name & sample == sample_name & rank <= n,
                    .(rank, species_label, rel_abundance_pct)]
  if (nrow(dt) == 0) return(character())
  paste0(
    "- ", method_name, " / ", sample_name, ": ",
    paste(sprintf("%s %.1f%%", dt$species_label, dt$rel_abundance_pct), collapse = "; ")
  )
}

samples <- sort(unique(results$sample))
report <- c(
  "# minitax bacteria RefSeq pilot analysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Output files",
  "",
  "- `minitax_bacteria_refseq_5_BestAln_taxa.all.sum.tsv`",
  "- `minitax_bacteria_refseq_5_SpeciesEstimate_taxa.all.sum.tsv`",
  "- `analysis/sample_totals_by_method.tsv`",
  "- `analysis/top20_species_by_method.tsv`",
  "- `analysis/top20_genera_by_method.tsv`",
  "- `analysis/species_method_relative_abundance_comparison.tsv`",
  "",
  "## Sample totals",
  "",
  capture.output(print(sample_totals)),
  "",
  "## Dominant species",
  ""
)
for (sample_id in samples) {
  report <- c(report, paste0("### ", sample_id), fmt_top("BestAln", sample_id), fmt_top("SpeciesEstimate", sample_id), "")
}

writeLines(report, file.path(analysis_dir, "minitax_5_result_summary.md"))

message("Wrote TSV exports and analysis files to: ", normalizePath(analysis_dir))
