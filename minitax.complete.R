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

message("Welcome to", "\t", "Ă„â€ÄąÄ„Ă˘â‚¬ÂÄąË‡Ă„ĹąĂ„ËťĂ‚Â©Ă„ĹąĂ„ËťĂ‚Â®Ă„ĹąĂ„ËťĂ‚Â©Ă„ĹąĂ„ËťĂ‚Â´Ä‚Ĺ˝Ă˘â‚¬ĹźÄ‚Ĺ˝Ă‚Â§", " !", "\n", "minitax, version: ", version)


message(
  "Ä‚ËĂ˘â‚¬ËĂ˘â‚¬ĹĄÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ˘â‚¬â€ť
Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â                                 Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â
Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â  Ä‚ËĂ˘â‚¬ĹĄÄąĹˇÄ‚ËĂ˘â‚¬ĹĄĂ‚Â¬Ä‚ËĂ˘â‚¬ĹĄĂ‚Â  Ä‚ËĂ˘â‚¬ĹĄĂ‚Â¬  Ä‚ËĂ˘â‚¬ĹĄĂ‚Â Ä‚ËĂ˘â‚¬ĹĄÄąĹˇ  Ä‚ËĂ˘â‚¬ĹĄĂ‚Â¬  Ä‚ËĂ˘â‚¬ĹĄĂ˘â€šÂ¬Ä‚ËĂ˘â‚¬ĹĄĂ‚Â¬Ä‚ËĂ˘â‚¬ĹĄĂ˘â€šÂ¬  Ä‚ËĂ˘â‚¬ĹĄÄąĹˇÄ‚ËĂ˘â‚¬ĹĄĂ˘â€šÂ¬Ä‚ËĂ˘â‚¬ĹĄĂ‚Â  Ä‚ËĂ˘â‚¬ËĂ‹â€şÄ‚ËĂ˘â‚¬ËĂ‚Â±   Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â
Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â  Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬ĹˇÄ‚ËĂ˘â‚¬ĹĄĂ˘â‚¬ĹˇÄ‚ËĂ˘â‚¬ĹĄĂ˘â‚¬Ĺˇ  Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬Ĺˇ  Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬ĹˇÄ‚ËĂ˘â‚¬ĹĄÄąÄ„Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬Ĺˇ  Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬Ĺˇ   Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬Ĺˇ   Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬Ĺˇ Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬Ĺˇ  ||   Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â
Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â  Ä‚ËĂ˘â‚¬ĹĄĂ‚Â´ Ä‚ËĂ˘â‚¬ĹĄĂ‚Â´  Ä‚ËĂ˘â‚¬ĹĄĂ‚Â´  Ä‚ËĂ˘â‚¬ĹĄĂ‚Â Ä‚ËĂ˘â‚¬ĹĄĂ˘â‚¬ĹĄ  Ä‚ËĂ˘â‚¬ĹĄĂ‚Â´   Ä‚ËĂ˘â‚¬ĹĄĂ‚Â´   Ä‚ËĂ˘â‚¬ĹĄĂ‚Â´ Ä‚ËĂ˘â‚¬ĹĄĂ‚Â´  Ä‚ËĂ˘â‚¬ËĂ‚Â±Ä‚ËĂ˘â‚¬ËĂ‹â€ş   Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â
Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â                                 Ä‚ËĂ˘â‚¬ËĂ˘â‚¬Â
Ä‚ËĂ˘â‚¬ËÄąË‡Ä‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËĂ‚ÂÄ‚ËĂ˘â‚¬ËÄąÄ„
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
  '.*_best_alignments_w_taxa.tsv',
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

keep.max.cigar <- parse_logical_config(config$keep.max.cigar, default = TRUE)
CIGAR_points   <- separate(data.frame(stringr::str_split_1(config$CIGAR_points, ';\\s*')), col = 1, into = c('key', 'value'), sep = '\\s*=\\s*')
CIGAR_points$value <- as.integer(CIGAR_points$value)
if (any(is.na(CIGAR_points$value))) {
  stop("CIGAR_points contains non-integer values.", call. = FALSE)
}
CIGAR_points <- CIGAR_points %>% spread(key, value)


#### configs used
message('Configs used: ')
print(t(config))
##

####
message('Constructed metadata:')
print(metadata)

#### Set up future to use multiple cores  with future
if (Sys.info()[["sysname"]] == 'Windows') {
  plan(multisession, workers = nproc)
} else {
  plan(multicore, workers = nproc)
}

####
message('everything is ready! Starting minitax!')


#### RUN MINITAX ON ALIGNMENTS ####


#### 1. Best Alignment Method -->> best precision !
if (any(methods.to.use == 'BestAln')) {
  ###
  methods  <- 'BestAln'
  message('\n \n ', methods, ' method started, on: ', date())

  ### from tax identity results
  ### from tax identity results
  taxa_or_bamfile <- list.files(file.path(outdir,  'best_alignments_w_taxa'),
                                '.*_best_alignments_w_taxa.tsv', full.names = T)
  if(length(taxa_or_bamfile) > 0) {
    input <- 'taxa.DT'
  } else {
    ### from bafiles
    taxa_or_bamfile <- paste0(pardir, '/', bamfiles[])
    input <- 'bamfile'
  }

  message('Running minitax on the following files: \n',
          paste(taxa_or_bamfile, collapse = '\n'))

  plan(multisession)
  taxa.sum.BestAln <- rbindlist(
    future_lapply(taxa_or_bamfile[],
                  wrap.fun.complete,
                  config=config,
                  input=input,
                  #sample=gsub(pattern, '', gsub('.*/', '', bamfiles)),
                  methods=methods,
                  keep.max.cigar=keep.max.cigar, CIGAR_points=CIGAR_points,
                  mapq.filt=mapq.filt,

                  saveRDS=F,
                  steps=c('1.Import_Alns', '2.Refine', '3.Find_taxidentity', '4.Summarise'),
                  best.mapq=best.mapq,

                  db=db, db.data=db.data, db.uni.data=db.uni.data,
                  ranks = ranks # c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    ), fill=TRUE)
  taxa.sum <- taxa.sum.BestAln
  setDT(taxa.sum)

  ###
  saveRDS(taxa.sum, method_result_path(methods, '_taxa.all.sum.rds'))
  fwrite(taxa.sum, method_result_path(methods, '_taxa.all.sum.tsv'), sep = '\t', na = 'NA')

  message('output of ', methods, ': ')
  print(head(taxa.sum))

  #### Generate phyloseq objects ####
  if (HAS_PHYLOSEQ) {
  ps <- NULL
  try({
    sample_data        <- metadata
    sample_data$method <- methods

    setDF(taxa.sum)
    taxa.sum <- taxa.sum[,c('tax.identity', 'lineage', ranks,  'sample', 'count')]
    taxa.sum[is.na(taxa.sum$tax.identity), c('tax.identity', 'lineage', ranks)] <- 'unclassified'

    setDT(taxa.sum)
    otutab <- dcast.data.table(taxa.sum, lineage ~ sample, value.var = 'count', fill=0)
    otutab <- data.frame(otutab[,-1], row.names = otutab$lineage)

    setDF(taxa.sum)
    taxtab <- data.frame(unique(taxa.sum[,c('lineage', ranks, 'tax.identity')]))
    rownames(taxtab) <- taxtab$lineage

    taxtab <- taxtab[rownames(otutab), ]

    ps <- phyloseq(otu_table(as.matrix(otutab), taxa_are_rows = T),
                   tax_table(as.matrix(taxtab)),
                   sample_data(sample_data))

    sample_names(ps)   <- paste(sample_data(ps)$workflow, sample_data(ps)$db, sample_data(ps)$method, sample_names(ps), sep='_')

    message('Generated phyloseq object: ')
    print(ps)

    ###
    saveRDS(ps, method_result_path(methods, '_PS.rds'))
    ps.minitax.BA <- ps
  })
  if(is.null(ps)) stop("phyloseq object could not be created!")
  } else {
    message("Skipping phyloseq object for ", methods, " because phyloseq is not installed.")
  }

  message(methods, ' method finished, on: ', date(), '\n \n ')

}
####



#### 2. Random Alignment Method
if (any(methods.to.use == 'RandAln')) {
  ###
  methods  <- 'RandAln'
  message('\n \n ', methods, ' method started, on: ', date())

  ### from tax identity results
  taxa_or_bamfile <- list.files(file.path(outdir,  'best_alignments_w_taxa'),
                                '.*_best_alignments_w_taxa.tsv', full.names = T)
  if(length(taxa_or_bamfile) > 0) {
    input <- 'taxa.DT'
  } else {
    ### from bafiles
    taxa_or_bamfile <- paste0(pardir, '/', bamfiles[])
    input <- 'bamfile'
  }

  message('Running minitax on the following files: \n',
          paste(taxa_or_bamfile, collapse = '\n'))

  plan(multisession)
  taxa.sum.RandAln <- rbindlist(
    future_lapply(taxa_or_bamfile[],
                  wrap.fun.complete,
                  config=config,
                  input=input,
                  #sample=gsub(pattern, '', gsub('.*/', '', bamfiles)),
                  methods=methods,
                  keep.max.cigar=keep.max.cigar, CIGAR_points=CIGAR_points,
                  mapq.filt=mapq.filt,

                  saveRDS=F,
                  steps=c('1.Import_Alns', '2.Refine', '3.Find_taxidentity', '4.Summarise'),
                  best.mapq=best.mapq,

                  db=db, db.data=db.data, db.uni.data=db.uni.data,
                  ranks = ranks # c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    ), fill=TRUE)
  taxa.sum <- taxa.sum.RandAln
  setDT(taxa.sum)

  ###
  saveRDS(taxa.sum, method_result_path(methods, '_taxa.all.sum.rds'))
  fwrite(taxa.sum, method_result_path(methods, '_taxa.all.sum.tsv'), sep = '\t', na = 'NA')
  #taxa.sum <- readRDS(paste0(outdir, '/', outdir, '_', 'taxa.all.sum_', methods, '.rds'))

  message('output of ', methods, ': ')
  print(head(taxa.sum))

  #### Generate phyloseq objects ####
  if (HAS_PHYLOSEQ) {
   ps <- NULL
  try({
    sample_data        <- metadata
    sample_data$method <- methods

    setDF(taxa.sum)
    taxa.sum <- taxa.sum[,c('tax.identity', 'lineage', ranks,  'sample', 'count')]
    taxa.sum[is.na(taxa.sum$tax.identity), c('tax.identity', 'lineage', ranks)] <- 'unclassified'

    setDT(taxa.sum)
    otutab <- dcast.data.table(taxa.sum, lineage ~ sample, value.var = 'count', fill=0)
    otutab <- data.frame(otutab[,-1], row.names = otutab$lineage)

    setDF(taxa.sum)
    taxtab <- data.frame(unique(taxa.sum[,c('lineage', ranks, 'tax.identity')]))
    rownames(taxtab) <- taxtab$lineage

    taxtab <- taxtab[rownames(otutab), ]

    ps <- phyloseq(otu_table(as.matrix(otutab), taxa_are_rows = T),
                   tax_table(as.matrix(taxtab)),
                   sample_data(sample_data))

    sample_names(ps)   <- paste(sample_data(ps)$workflow, sample_data(ps)$db, sample_data(ps)$method, sample_names(ps), sep='_')

    message('Generated phyloseq object: ')
    print(ps)

    ###
    saveRDS(ps, method_result_path(methods, '_PS.rds'))
    ps.minitax.RA <- ps
  })
  if(is.null(ps)) stop("phyloseq object could not be created!")
  } else {
    message("Skipping phyloseq object for ", methods, " because phyloseq is not installed.")
  }

  message(methods, ' method finished, on: ', date(), '\n \n ')

}
####


#### 3. LCA Method  -->> Summarize the reads on their tax.identity
if (any(methods.to.use == 'LCA')) {
  ###
  methods  <- 'LCA'
  message(methods, ' method started, on: ', date())

  ### from tax identity results
  taxa_or_bamfile <- list.files(file.path(outdir,  'best_alignments_w_taxa'),
                                '.*_best_alignments_w_taxa.tsv', full.names = T)
  if(length(taxa_or_bamfile) > 0) {
    input <- 'taxa.DT'
  } else {
    ### from bafiles
    taxa_or_bamfile <- paste0(pardir, '/', bamfiles[])
    input <- 'bamfile'
  }

  message('Running minitax on the following files: \n',
          paste(taxa_or_bamfile, collapse = '\n'))

  plan(multisession)

  taxa.sum.LCA <- rbindlist(
     future_lapply(taxa_or_bamfile[],
                   wrap.fun.complete,
                   config=config,
                   input=input,
                   #sample=gsub(pattern, '', gsub('.*/', '', bamfiles)),
                   methods=methods,
                   keep.max.cigar=keep.max.cigar, CIGAR_points=CIGAR_points,
                   mapq.filt=mapq.filt,

                   saveRDS=F,
                   steps=c('1.Import_Alns', '2.Refine', '3.Find_taxidentity', '4.Summarise'),
                   best.mapq=best.mapq,

                   db=db, db.data=db.data, db.uni.data=db.uni.data,
                   ranks = ranks # c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    ), fill=TRUE)
  taxa.sum.LCA <- rename_duptaxa(taxa.sum.LCA, ranks=ranks)
  taxa.sum <- taxa.sum.LCA


  ###
  saveRDS(taxa.sum, method_result_path(methods, '_taxa.all.sum.rds'))
  fwrite(taxa.sum, method_result_path(methods, '_taxa.all.sum.tsv'), sep = '\t', na = 'NA')
  #taxa.sum <- readRDS(paste0(outdir, '/', outdir, '_', 'taxa.all.sum_', methods, '.rds'))

  message('output of ', methods, ': ')
  print(head(taxa.sum))


  #### Generate phyloseq objects ####
  if (HAS_PHYLOSEQ) {
  ps <- NULL
  try({
    sample_data        <- metadata
    sample_data$method <- methods

    setDF(taxa.sum)
    taxa.sum <- taxa.sum[,c('tax.identity', 'lineage', ranks,  'sample', 'count')]
    taxa.sum[is.na(taxa.sum$tax.identity), c('tax.identity', 'lineage', ranks)] <- 'unclassified'

    setDT(taxa.sum)
    otutab <- dcast.data.table(taxa.sum, lineage ~ sample, value.var = 'count', fill=0)
    otutab <- data.frame(otutab[,-1], row.names = otutab$lineage)

    setDF(taxa.sum)
    taxtab <- data.frame(unique(taxa.sum[,c('lineage', ranks, 'tax.identity')]))
    rownames(taxtab) <- taxtab$lineage

  taxtab <- taxtab[rownames(otutab), ]

    ps <- phyloseq(otu_table(as.matrix(otutab), taxa_are_rows = T),
                   tax_table(as.matrix(taxtab)),
                   sample_data(sample_data))

    sample_names(ps)   <- paste(sample_data(ps)$workflow, sample_data(ps)$db, sample_data(ps)$method, sample_names(ps), sep='_')

    message('Generated phyloseq object: ')
    print(ps)

    ###
    saveRDS(ps, method_result_path(methods, '_PS.rds'))
    ps.minitax.LCA <- ps
  })
  if(is.null(ps)) stop("phyloseq object could not be created!")
  } else {
    message("Skipping phyloseq object for ", methods, " because phyloseq is not installed.")
  }

  message(methods, ' method finished, on: ', date(), '\n \n ')

}
####


#### 4. Species Estimation Method -->> much faster !
if (any(methods.to.use == 'SpeciesEstimate')) {
  ###
  methods  <- 'SpeciesEstimate'
  message(methods, ' method started, on: ', date())

  ### from tax identity results
  taxa_or_bamfile <- list.files(file.path(outdir,  'best_alignments_w_taxa'),
                                '.*_best_alignments_w_taxa.tsv', full.names = T)
  if(length(taxa_or_bamfile) > 0) {
    input <- 'taxa.DT'
  } else {
    ### from bafiles
    taxa_or_bamfile <- paste0(pardir, '/', bamfiles[])
    input <- 'bamfile'
  }

  message('Running minitax on the following files: \n',
          paste(taxa_or_bamfile, collapse = '\n'))

  plan(multisession)

  taxID <- ranks[length(ranks)]

  taxa.sum.SpecEst <- rbindlist(
    future_lapply(taxa_or_bamfile[],
                  wrap.fun.complete,
                  config=config,
                  input=input,
                  #sample=gsub(pattern, '', gsub('.*/', '', bamfiles)),
                  methods=methods,
                  keep.max.cigar=keep.max.cigar, CIGAR_points=CIGAR_points,
                  mapq.filt=mapq.filt,

                  saveRDS=F,
                  steps=c('1.Import_Alns', '2.Refine', '3.Find_taxidentity', '4.Summarise'),
                  best.mapq=best.mapq,

                  db=db, db.data=db.data, db.uni.data=db.uni.data,
                  ranks = ranks # c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    ), fill=TRUE)

  taxa.sum.SpecEst[, tax.identity := get(taxID)]
  taxa.sum.SpecEst[is.na(tax.identity), tax.identity:= 'unclassified']
  taxa.sum.SpecEst[, lineage  := tax.identity]

  taxa.sum <- taxa.sum.SpecEst

  ###
  saveRDS(taxa.sum, method_result_path(methods, '_taxa.all.sum.rds'))
  fwrite(taxa.sum, method_result_path(methods, '_taxa.all.sum.tsv'), sep = '\t', na = 'NA')
  #taxa.sum <- readRDS(paste0(outdir, '/', outdir, '_', 'taxa.all.sum_', methods, '.rds'))

  message('output of ', methods, ': ')
  print(head(taxa.sum))



  #### Generate phyloseq objects ####
  if (HAS_PHYLOSEQ) {
  ps <- NULL
  try({

    ## for species estimate, the output is normalized
    taxa.sum[,count := norm_count]

    #
    sample_data        <- metadata
    sample_data$method <- methods

    setDF(taxa.sum)
    taxa.sum <- taxa.sum[,c('tax.identity', 'lineage', ranks,  'sample', 'count')]
    taxa.sum[is.na(taxa.sum$tax.identity), c('tax.identity', 'lineage', ranks)] <- 'unclassified'

    setDT(taxa.sum)
    otutab <- dcast.data.table(taxa.sum, lineage ~ sample, value.var = 'count', fill=0)
    otutab <- data.frame(otutab[,-1], row.names = otutab$lineage)

    setDF(taxa.sum)
    taxtab <- data.frame(unique(taxa.sum[,c('lineage', ranks, 'tax.identity')]))
    rownames(taxtab) <- taxtab$lineage

    taxtab <- taxtab[rownames(otutab), ]

    ps <- phyloseq(otu_table(as.matrix(otutab), taxa_are_rows = T),
                   tax_table(as.matrix(taxtab)),
                   sample_data(sample_data))

    sample_names(ps)   <- paste(sample_data(ps)$workflow, sample_data(ps)$db, sample_data(ps)$method, sample_names(ps), sep='_')

    message('Generated phyloseq object: ')
    print(ps)

    ###
    saveRDS(ps, method_result_path(methods, '_PS.rds'))
    ps.minitax.SE <- ps
  })
  if(is.null(ps)) stop("phyloseq object could not be created!")
  } else {
    message("Skipping phyloseq object for ", methods, " because phyloseq is not installed.")
  }

  message(methods, ' method finished, on: ', date(), '\n \n ')

}
####
#### ####

##############################
