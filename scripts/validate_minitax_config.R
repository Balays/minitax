#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Validate a minitax config before mapping or classification.\n\n",
    "Usage:\n",
    "  scripts/validate_minitax_config.R CONFIG [--step all|map|classify|database]\n\n",
    "Options:\n",
    "  --step STEP   Check only one workflow part. Default: all\n",
    "  -h, --help    Show this help.\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = ifelse(length(args) == 0, 1, 0))
}

config_file <- NULL
step <- "all"
i <- 1
while (i <= length(args)) {
  arg <- args[[i]]
  if (arg == "--step") {
    if (i == length(args)) stop("--step requires a value.", call. = FALSE)
    step <- args[[i + 1]]
    i <- i + 2
  } else if (startsWith(arg, "--")) {
    stop("Unknown option: ", arg, call. = FALSE)
  } else if (is.null(config_file)) {
    config_file <- arg
    i <- i + 1
  } else {
    stop("Unexpected argument: ", arg, call. = FALSE)
  }
}

allowed_steps <- c("all", "map", "classify", "database")
if (!step %in% allowed_steps) {
  stop("--step must be one of: ", paste(allowed_steps, collapse = ", "), call. = FALSE)
}
if (is.null(config_file)) {
  usage()
  stop("Missing CONFIG.", call. = FALSE)
}

checks <- data.frame(level = character(), message = character(), stringsAsFactors = FALSE)
`%+%` <- function(a, b) paste0(a, b)

add_check <- function(level, message) {
  checks <<- rbind(checks, data.frame(level = level, message = message, stringsAsFactors = FALSE))
  cat(sprintf("[%s] %s\n", level, message))
}
ok <- function(message) add_check("OK", message)
warn <- function(message) add_check("WARN", message)
err <- function(message) add_check("ERROR", message)
info <- function(message) add_check("INFO", message)

is_config_na <- function(value) {
  length(value) == 0 || is.na(value) || !nzchar(trimws(value)) || toupper(trimws(value)) %in% c("NA", "NULL")
}

trim_quotes <- function(value) {
  value <- trimws(value)
  sub("^(['\"])(.*)\\1$", "\\2", value)
}

has_surrounding_quotes <- function(value) {
  grepl("^(['\"]).*\\1$", trimws(value))
}

split_values <- function(value, default = character()) {
  if (is_config_na(value)) return(default)
  values <- trimws(unlist(strsplit(value, "[;,]")))
  values[nzchar(values)]
}

normalize_config_path <- function(path) {
  if (is_config_na(path)) return(NA_character_)
  path <- trim_quotes(path)
  sys <- Sys.info()[["sysname"]]
  if (!is.na(sys) && sys == "Windows") {
    path <- gsub("^/mnt/([A-Za-z])/", "\\U\\1:/", path, perl = TRUE)
  } else {
    path <- gsub("^([A-Za-z]):[\\\\/]", "/mnt/\\L\\1/", path, perl = TRUE)
    path <- gsub("\\\\", "/", path)
  }
  path
}

script_path <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
script_dir <- if (!is.na(script_path) && nzchar(script_path)) {
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  file.path(getwd(), "scripts")
}
repo_dir <- dirname(script_dir)

read_config <- function(path) {
  if (!file.exists(path)) stop("Config file not found: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  lines <- sub("\r$", "", lines)
  config <- character()
  duplicates <- character()
  for (line in lines) {
    if (!nzchar(trimws(line)) || startsWith(trimws(line), "#")) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 2) next
    key <- trimws(fields[[1]])
    value <- trimws(fields[[2]])
    if (!nzchar(key) || key == "argument") next
    if (key %in% names(config)) duplicates <- c(duplicates, key)
    config[[key]] <- value
  }
  attr(config, "duplicates") <- unique(duplicates)
  config
}

cfg <- read_config(config_file)
info("Config file: " %+% normalizePath(config_file, mustWork = FALSE))
info("Validation step: " %+% step)
info("Working directory: " %+% getwd())
if (length(attr(cfg, "duplicates")) > 0) {
  warn("Duplicate config key(s); last value wins: " %+% paste(attr(cfg, "duplicates"), collapse = ", "))
}

get_cfg <- function(name, default = NA_character_) {
  if (!name %in% names(cfg) || is_config_na(cfg[[name]])) return(default)
  cfg[[name]]
}

require_keys <- function(keys) {
  missing_by_name <- keys[!keys %in% names(cfg)]
  present <- intersect(keys, names(cfg))
  missing_by_value <- present[vapply(present, function(k) is_config_na(cfg[[k]]), logical(1))]
  missing <- unique(c(missing_by_name, missing_by_value))
  if (length(missing) > 0) {
    err("Missing required config argument(s): " %+% paste(missing, collapse = ", "))
  } else {
    ok("Required config arguments are present: " %+% paste(keys, collapse = ", "))
  }
}

check_no_quotes <- function(name) {
  value <- get_cfg(name)
  if (!is_config_na(value) && has_surrounding_quotes(value)) {
    err("Config value '" %+% name %+% "' has surrounding quotes. Use " %+% trim_quotes(value) %+% " instead.")
  }
}

check_choice <- function(name, choices, required = TRUE) {
  value <- get_cfg(name)
  if (is_config_na(value)) {
    if (required) err("Missing config value: " %+% name)
    return(invisible(FALSE))
  }
  if (!value %in% choices) {
    err("Config value '" %+% name %+% "' must be one of: " %+% paste(choices, collapse = ", "))
    return(invisible(FALSE))
  }
  ok(name %+% " = " %+% value)
  invisible(TRUE)
}

check_integer <- function(name, min_value = NULL, max_value = NULL, required = TRUE) {
  value <- get_cfg(name)
  if (is_config_na(value)) {
    if (required) err("Missing integer config value: " %+% name)
    return(invisible(NA_integer_))
  }
  number <- suppressWarnings(as.integer(value))
  if (is.na(number) || as.character(number) != trimws(value)) {
    err(name %+% " must be an integer.")
    return(invisible(NA_integer_))
  }
  if (!is.null(min_value) && number < min_value) err(name %+% " must be >= " %+% min_value %+% ".")
  if (!is.null(max_value) && number > max_value) err(name %+% " must be <= " %+% max_value %+% ".")
  if ((is.null(min_value) || number >= min_value) && (is.null(max_value) || number <= max_value)) {
    ok(name %+% " = " %+% number)
  }
  invisible(number)
}

check_index_batch_size <- function(name, default = "128G") {
  value <- get_cfg(name, default)
  if (is_config_na(value)) value <- default
  if (!grepl("^[0-9]+[KkMmGgTt]?$", value)) {
    err(name %+% " must look like 128G, 64000M, or 500000000.")
  } else {
    ok(name %+% " = " %+% value)
  }
}

check_logical <- function(name, required = FALSE) {
  value <- get_cfg(name)
  if (is_config_na(value)) {
    if (required) err("Missing logical config value: " %+% name)
    return(invisible(FALSE))
  }
  valid <- toupper(trimws(value)) %in% c("TRUE", "T", "YES", "Y", "1", "FALSE", "F", "NO", "N", "0")
  if (!valid) err(name %+% " must be TRUE/FALSE, T/F, yes/no, or 1/0.") else ok(name %+% " is a valid logical value.")
  invisible(valid)
}

path_exists <- function(path, label, required = TRUE) {
  path <- normalize_config_path(path)
  if (is_config_na(path)) {
    if (required) err(label %+% " is missing.")
    return(invisible(FALSE))
  }
  if (file.exists(path)) {
    ok(label %+% " exists: " %+% path)
    return(invisible(TRUE))
  }
  if (required) err(label %+% " not found: " %+% path) else warn(label %+% " not found: " %+% path)
  invisible(FALSE)
}

file_readable <- function(path, label, required = TRUE) {
  path <- normalize_config_path(path)
  if (!path_exists(path, label, required = required)) return(invisible(FALSE))
  if (file.access(path, 4) == 0) {
    ok(label %+% " is readable.")
    invisible(TRUE)
  } else {
    if (required) err(label %+% " is not readable: " %+% path) else warn(label %+% " is not readable: " %+% path)
    invisible(FALSE)
  }
}

check_executable <- function(path, label, required = TRUE) {
  path <- normalize_config_path(path)
  if (!is_config_na(path) && !grepl("[/\\\\]", path)) {
    return(check_command(path, label, required = required))
  }
  if (!path_exists(path, label, required = required)) return(invisible(FALSE))
  if (file.access(path, 1) == 0) {
    ok(label %+% " is executable.")
    invisible(TRUE)
  } else {
    warn(label %+% " exists but is not marked executable: " %+% path)
    invisible(TRUE)
  }
}

check_command <- function(command, label, required = TRUE) {
  found <- Sys.which(command)
  if (nzchar(found)) {
    ok(label %+% " found: " %+% found)
    invisible(TRUE)
  } else {
    if (required) err(label %+% " not found in PATH.") else warn(label %+% " not found in PATH.")
    invisible(FALSE)
  }
}

check_package <- function(pkg, optional = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    ok("R package installed: " %+% pkg)
  } else if (optional) {
    warn("Optional R package not installed: " %+% pkg)
  } else {
    err("Required R package not installed: " %+% pkg)
  }
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

list_fastq <- function(indir, suffix) {
  if (!dir.exists(indir)) return(character())
  files <- list.files(indir, full.names = FALSE, no.. = TRUE)
  files[endsWith(files, suffix)]
}

check_fastq_inputs <- function() {
  indir <- normalize_config_path(get_cfg("indir"))
  suffix <- get_cfg("fastq_suffix")
  platform <- get_cfg("platform")
  reads <- get_cfg("reads", "paired")

  path_exists(indir, "FASTQ input directory", required = TRUE)
  if (is_config_na(suffix)) {
    err("fastq_suffix is required for mapping.")
    return(invisible(FALSE))
  }
  if (!dir.exists(indir)) return(invisible(FALSE))

  fastqs <- list_fastq(indir, suffix)
  if (length(fastqs) == 0) {
    err("No FASTQ files ending with '" %+% suffix %+% "' found in " %+% indir)
    return(invisible(FALSE))
  }

  if (platform == "Illumina" && reads == "paired") {
    pattern <- get_cfg("fastq_pair_pattern", "")
    if (is_config_na(pattern)) pattern <- ""
    r1_suffix <- pattern %+% "_R1" %+% suffix
    r2_suffix <- pattern %+% "_R2" %+% suffix
    r1 <- fastqs[endsWith(fastqs, r1_suffix)]
    r2 <- fastqs[endsWith(fastqs, r2_suffix)]
    if (length(r1) == 0) err("No Illumina R1 FASTQ files found with suffix " %+% r1_suffix)
    if (length(r2) == 0) err("No Illumina R2 FASTQ files found with suffix " %+% r2_suffix)
    bases_r1 <- sub(escape_regex(r1_suffix) %+% "$", "", r1)
    bases_r2 <- sub(escape_regex(r2_suffix) %+% "$", "", r2)
    missing_r2 <- setdiff(bases_r1, bases_r2)
    missing_r1 <- setdiff(bases_r2, bases_r1)
    if (length(missing_r2) > 0) err("Missing R2 pair(s) for sample(s): " %+% paste(missing_r2, collapse = ", "))
    if (length(missing_r1) > 0) err("Missing R1 pair(s) for sample(s): " %+% paste(missing_r1, collapse = ", "))
    if (length(r1) > 0 && length(r2) > 0 && length(missing_r1) == 0 && length(missing_r2) == 0) {
      ok("Illumina paired FASTQ files look matched for " %+% length(r1) %+% " sample(s).")
    }
  } else {
    ok("Found " %+% length(fastqs) %+% " FASTQ file(s) ending with " %+% suffix)
  }
}

check_mapq <- function() {
  value <- get_cfg("mapq.filt")
  if (is_config_na(value)) {
    ok("mapq.filt = NA")
    return(invisible(TRUE))
  }
  value <- trimws(value)
  if (!grepl("^[0-9]+(:[0-9]+)?$", value)) {
    err("mapq.filt must be NA, an integer, or an integer range like 0:59.")
    return(invisible(FALSE))
  }
  parts <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
  if (length(parts) == 2 && parts[[1]] > parts[[2]]) {
    err("mapq.filt range start must be <= range end.")
  } else {
    ok("mapq.filt is valid: " %+% value)
  }
}

check_cigar_points <- function() {
  value <- get_cfg("CIGAR_points")
  if (is_config_na(value)) {
    warn("CIGAR_points is missing; minitax will use internal defaults.")
    return(invisible(TRUE))
  }
  pieces <- trimws(unlist(strsplit(value, ";", fixed = TRUE)))
  pieces <- pieces[nzchar(pieces)]
  parsed <- strsplit(pieces, "\\s*=\\s*")
  keys <- vapply(parsed, `[`, character(1), 1)
  vals <- suppressWarnings(as.integer(vapply(parsed, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))))
  required <- c("match_score", "mismatch_score", "insertion_score", "deletion_score", "gap_opening_penalty", "gap_extension_penalty")
  missing <- setdiff(required, keys)
  if (length(missing) > 0) err("CIGAR_points missing key(s): " %+% paste(missing, collapse = ", "))
  if (any(is.na(vals))) err("CIGAR_points values must be integers.")
  if (length(missing) == 0 && !any(is.na(vals))) ok("CIGAR_points has all required integer scores.")
}

check_methods <- function() {
  methods <- split_values(get_cfg("methods", "BestAln"), default = "BestAln")
  allowed <- c("BestAln", "RandAln", "LCA", "SpeciesEstimate")
  bad <- setdiff(methods, allowed)
  if (length(bad) > 0) err("Unknown method(s): " %+% paste(bad, collapse = ", "))
  if (length(bad) == 0) ok("Method list is valid: " %+% paste(methods, collapse = ", "))

  thresholds <- split_values(get_cfg("BestAln_thresholds"), default = "0.6")
  nums <- suppressWarnings(as.numeric(thresholds))
  if (any(is.na(nums))) {
    err("BestAln_thresholds must be numeric.")
  } else if (!length(nums) %in% c(1, 7)) {
    err("BestAln_thresholds must contain either one value or seven rank-level values.")
  } else if (any(nums < 0 | nums > 1)) {
    warn("BestAln_thresholds are usually probabilities between 0 and 1.")
  } else {
    ok("BestAln_thresholds is valid.")
  }
}

check_output_dir <- function() {
  db <- get_cfg("db")
  platform <- get_cfg("platform", "platform")
  project <- get_cfg("project", "project")
  vregion <- get_cfg("Vregion", "WGS")
  outdir <- get_cfg("outdir")
  if (is_config_na(outdir)) {
    outdir <- paste("minitax", project, platform, vregion, db, sep = "_")
  }
  outdir <- normalize_config_path(outdir)
  parent <- dirname(outdir)
  if (!dir.exists(parent)) {
    warn("Output parent directory does not exist yet: " %+% parent)
  } else if (file.access(parent, 2) != 0) {
    err("Output parent directory is not writable: " %+% parent)
  } else {
    ok("Output parent directory is writable: " %+% parent)
  }
  if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)) > 0) {
    warn("Output directory already exists and is not empty: " %+% outdir)
  } else {
    ok("Output directory is available: " %+% outdir)
  }
}

check_disk_space <- function() {
  outdir <- get_cfg("outdir")
  if (is_config_na(outdir)) outdir <- "."
  outdir <- normalize_config_path(outdir)
  probe <- if (dir.exists(outdir)) outdir else dirname(outdir)
  if (!dir.exists(probe) || !nzchar(Sys.which("df"))) return(invisible(FALSE))
  result <- tryCatch(system2("df", c("-Pk", probe), stdout = TRUE, stderr = FALSE), error = function(e) character())
  if (length(result) < 2) return(invisible(FALSE))
  fields <- strsplit(gsub("^\\s+|\\s+$", "", result[[2]]), "\\s+")[[1]]
  if (length(fields) >= 4) {
    free_gb <- suppressWarnings(as.numeric(fields[[4]]) / 1024 / 1024)
    if (!is.na(free_gb)) {
      if (free_gb < 10) warn(sprintf("Low free disk space near output path: %.1f GB", free_gb)) else ok(sprintf("Free disk space near output path: %.1f GB", free_gb))
    }
  }
}

check_database <- function() {
  db <- get_cfg("db")
  db_dir <- normalize_config_path(get_cfg("db.dir"))
  file_readable(db_dir, "Database directory", required = TRUE)
  if (!dir.exists(db_dir)) return(invisible(FALSE))

  if (db == "all_NCBI_genomes") {
    required_files <- c("NCBI.db.tsv", "NCBI.db.uni.tsv", "NCBI.db.genomesize.tsv")
    optional_files <- c("NCBI.db.uni.spec.tsv", "NCBI_genome_collection_seqlengths.txt")
  } else if (db == "EMUdb") {
    required_files <- c("taxonomy.tsv", "species_taxid.fasta")
    optional_files <- character()
  } else {
    required_files <- c("db_data.tsv")
    optional_files <- character()
    warn("Unknown/custom database type '" %+% db %+% "'. Expecting generic db_data.tsv.")
  }

  for (file in required_files) {
    file_readable(file.path(db_dir, file), "Database file " %+% file, required = TRUE)
  }
  for (file in optional_files) {
    file_readable(file.path(db_dir, file), "Optional database file " %+% file, required = FALSE)
  }

  index <- get_cfg("mm2_index")
  if (!is_config_na(index)) {
    index_path <- file.path(db_dir, index)
    index_required <- step %in% c("all", "map")
    if (file.exists(index_path)) {
      file_readable(index_path, "minimap2 index", required = index_required)
    } else {
      mm2_ref <- normalize_config_path(get_cfg("mm2_ref"))
      ref_candidates <- character()
      if (!is_config_na(mm2_ref)) {
        ref_candidates <- unique(c(mm2_ref, file.path(db_dir, mm2_ref)))
      }
      ref_found <- ref_candidates[file.exists(ref_candidates)]
      if (index_required && length(ref_found) > 0) {
        warn("minimap2 index not found, but mm2_ref exists; minitax.sh can build the index: " %+% ref_found[[1]])
      } else {
        file_readable(index_path, "minimap2 index", required = index_required)
      }
    }
  } else if (step %in% c("all", "map")) {
    err("mm2_index is required for mapping.")
  }
}

check_r_sources <- function() {
  minitax_dir <- get_cfg("minitax.dir")
  misc_dir <- get_cfg("misc.dir")
  if (is_config_na(minitax_dir)) minitax_dir <- repo_dir
  if (is_config_na(misc_dir)) misc_dir <- repo_dir
  minitax_dir <- normalize_config_path(minitax_dir)
  misc_dir <- normalize_config_path(misc_dir)

  required_sources <- file.path(minitax_dir, c(
    "R/minitax.v2.R",
    "R/cigar_score.R",
    "R/cigar_to_length.R",
    "R/cigar.sum.R",
    "R/rename_duptaxa.R",
    "R/minitax.wrapfun.complete.R",
    "R/BestAln.R",
    "R/add_lineage.R"
  ))
  misc_sources <- file.path(misc_dir, c(
    "R/ov.from.bam2.R",
    "R/dt.from.bam.R",
    "R/get.best.aln.R",
    "R/misc.functions.R",
    "R/add_rightmost_non.NA_col.R",
    "R/paste_columns_dt.R",
    "R/bam.flags.tsv"
  ))
  missing <- c(required_sources, misc_sources)[!file.exists(c(required_sources, misc_sources))]
  if (length(missing) == 0) {
    ok("All required R source/helper files are present.")
  } else {
    err("Missing R source/helper file(s): " %+% paste(missing, collapse = ", "))
  }
}

check_classify_inputs <- function() {
  db <- get_cfg("db")
  required_packages <- c("Rsamtools", "readr", "ggplot2", "tidyr", "dplyr", "GenomicAlignments", "data.table", "stringi", "stringr", "future.apply")
  if (db == "EMUdb") required_packages <- c(required_packages, "seqinr")
  for (pkg in required_packages) check_package(pkg)
  check_package("phyloseq", optional = TRUE)
  check_r_sources()

  if (step == "classify") {
    outdir <- get_cfg("outdir")
    if (is_config_na(outdir)) {
      outdir <- paste("minitax", get_cfg("project", "project"), get_cfg("platform", "platform"), get_cfg("Vregion", "WGS"), get_cfg("db", "db"), sep = "_")
    }
    outdir <- normalize_config_path(outdir)
    bam_dir <- file.path(outdir, "bam")
    taxa_dir <- file.path(outdir, "best_alignments_w_taxa")
    bam_files <- if (dir.exists(bam_dir)) list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE) else character()
    taxa_files <- if (dir.exists(taxa_dir)) list.files(taxa_dir, pattern = "_best_alignments_w_taxa\\.tsv$", full.names = TRUE) else character()
    if (length(bam_files) == 0 && length(taxa_files) == 0) {
      err("No BAM files or cached best_alignments_w_taxa files found for classification under: " %+% outdir)
    } else {
      ok("Classification inputs found: " %+% length(bam_files) %+% " BAM file(s), " %+% length(taxa_files) %+% " cached taxa file(s).")
    }
  }
}

require_keys(c("platform", "db", "db.dir", "nproc"))
for (name in c("db", "db.dir", "indir", "outdir", "mm2_path", "mm2_index", "mm2_ref", "mapper_backend", "parabricks_image", "mm2_index_batch", "minitax.dir", "misc.dir")) {
  check_no_quotes(name)
}
check_choice("platform", c("Illumina", "ONT", "PacBio"))
mapper_backend <- get_cfg("mapper_backend", "auto")
if (!mapper_backend %in% c("auto", "mm2-fast", "minimap2", "cpu", "parabricks")) {
  err("mapper_backend must be one of: auto, mm2-fast, minimap2, cpu, parabricks.")
} else {
  ok("mapper_backend = " %+% mapper_backend)
}
check_index_batch_size("mm2_index_batch")
parabricks_num_gpus <- get_cfg("parabricks_num_gpus", "1")
if (!grepl("^[0-9]+$", parabricks_num_gpus) || as.integer(parabricks_num_gpus) < 1) {
  err("parabricks_num_gpus must be a positive integer.")
} else {
  ok("parabricks_num_gpus = " %+% parabricks_num_gpus)
}
check_integer("nproc", min_value = 1)
check_integer("dt.threads", min_value = 1, required = FALSE)
cores <- parallel::detectCores(logical = TRUE)
if (!is.na(cores)) {
  nproc <- suppressWarnings(as.integer(get_cfg("nproc")))
  if (!is.na(nproc) && nproc > cores) warn("nproc (" %+% nproc %+% ") is greater than detected cores (" %+% cores %+% ").")
}
for (name in c("debug", "keep.highest.mapq.aln.only", "crop.na.tax", "multicore", "saveRAM", "reuse.taxa.cache", "keep.max.cigar", "best.mapq")) {
  check_logical(name, required = FALSE)
}
check_mapq()
check_cigar_points()
check_methods()
check_output_dir()
check_disk_space()

if (step %in% c("all", "map", "database", "classify")) {
  check_database()
}

if (step %in% c("all", "map")) {
  mm2_path <- get_cfg("mm2_path")
  mapper_backend <- get_cfg("mapper_backend", "auto")
  if (is_config_na(mm2_path) && mapper_backend != "parabricks") {
    err("mm2_path is required by minitax.sh unless mapper_backend=parabricks.")
  } else if (is_config_na(mm2_path)) {
    warn("mm2_path is missing; mapper_backend=parabricks will not be able to build a missing minimap2 index.")
  } else {
    check_executable(mm2_path, "minimap2 executable", required = TRUE)
  }
  if (mapper_backend == "parabricks") {
    mm2_ref <- get_cfg("mm2_ref")
    if (is_config_na(mm2_ref)) {
      err("mapper_backend=parabricks requires mm2_ref to point to the reference FASTA.")
    }
    check_command("docker", "docker", required = FALSE)
    check_command("nvidia-smi", "nvidia-smi", required = FALSE)
  }
  check_command("samtools", "samtools", required = TRUE)
  check_integer("Nsec", min_value = 0)
  platform <- get_cfg("platform")
  if (platform == "Illumina") {
    check_choice("reads", c("paired", "merged"))
  }
  check_fastq_inputs()
}

if (step %in% c("all", "classify")) {
  check_classify_inputs()
}

errors <- sum(checks$level == "ERROR")
warnings <- sum(checks$level == "WARN")
cat("\nSummary: ", errors, " error(s), ", warnings, " warning(s).\n", sep = "")
if (errors > 0) {
  cat("Config validation failed.\n")
  quit(status = 1)
}
cat("Config validation passed.\n")
