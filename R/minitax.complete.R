#!/usr/bin/env Rscript

#### Config ####
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=0) {
  config <- args
}
#### ####
## 

#### packages #####
library(Rsamtools, quietly = T)
library(readr, quietly = T)
library(ggplot2, quietly = T)
library(tidyr, quietly = T)
library(dplyr, quietly = T)
library(GenomicAlignments, quietly = T)
library(data.table, quietly = T)
library(stringi, quietly = T)
library(stringr, quietly = T)
library(phyloseq, quietly = T)
library(future.apply, quietly = T)
library(data.table, quietly = T)

#### ####


#### ####
config <- read.delim(config)
#print(config[,c(1:2)])

config <- config[,c(1:2)] %>% spread(argument, value)
t(config)

nproc <- as.integer(config$nproc) #args[2] #32 # 
#### ####
#stop()

#### OS #####
if (data.frame(Sys.info())[1,1] == 'Windows') {
  config$minitax.dir  <- gsub('/mnt/e/', 'E:/', config$minitax.dir)
  config$misc.dir     <- gsub('/mnt/e/', 'E:/', config$misc.dir)
  config$db.dir       <- gsub('/mnt/e/', 'E:/', config$db.dir)
}
#### ####
## 

#### functions ####
## minitax
minitax.dir <- config$minitax.dir ## 'my.R.packages/minitax'
source(paste0( minitax.dir, '/R/minitax.v2.R'))
source(paste0( minitax.dir, '/R/cigar_score.R'))
source(paste0( minitax.dir, '/R/cigar_to_length.R'))
source(paste0( minitax.dir, '/R/PS_from_taxa.sum.R'))
source(paste0( minitax.dir, '/R/minitax.R'))
source(paste0( minitax.dir, '/R/tax.identity.ranks.R'))
source(paste0( minitax.dir, '/R/make.krona.out.R'))
source(paste0( minitax.dir, '/R/estimate.species.R'))
source(paste0( minitax.dir, '/R/import.bam.R'))
source(paste0( minitax.dir, '/R/cigar.sum.R'))
source(paste0( minitax.dir, '/R/get_chunks.R'))
source(paste0( minitax.dir, '/R/taxprofile.from.minitax.R'))
source(paste0( minitax.dir, '/R/rename_duptaxa.R'))
source(paste0( minitax.dir, '/R/minitax.wrapfun.complete.R'))


## misc package
misc.dir <- config$misc.dir ## 'my.R.packages/Rlyeh-main'
source(paste0( misc.dir, '/ov.from.bam2.R'))
source(paste0( misc.dir, '/get.best.aln.R'))
source(paste0( misc.dir, '/misc.functions.R'))
bam.flags <- read.delim(paste0( misc.dir, '/bam.flags.tsv'))
#### ####
##


#### Project options ####
db        <- config$db ##'proGcontigs'
db.dir    <- config$db.dir ## paste0(path.prefix, 'data/databases/proGenomes')
project   <- config$project ## 'MCM'
Vregion   <- config$Vregion ## 'v1_2'
platform  <- config$platform ##'illumina'
outdir    <- paste('minitax', project, platform, Vregion, db, sep = '_')
#outdir <- paste0(outdir, '.partII')
dir.create(outdir)


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
}



message(db.name, ' database loaded')
#### ####
##


#### Metadata and parent directory of bamfiles from minimap2 output ####
pardir    <- paste0(outdir, '/bam')
pattern   <- '.bam'
bamfiles  <- grep('.bai', list.files(pardir, pattern = pattern), value = T, invert = T)


## filter
#bamfiles <- bamfiles[14]

samples   <- gsub(pattern, '', gsub('.*/', '', bamfiles))

#source('make.metadata.WGS.R')
#metadata$db <- db

metadata <- data.frame(sample=samples, bamfile=bamfiles, 
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

message('everything is ready!')


##### minitax v2 #####

#bamfile <- paste0(pardir, '/', sample, pattern)
#message('running minitax on ', bamfile, '...')

#### Minitax settings
keep.highest.mapq.aln.only <- config$keep.highest.mapq.aln.only
crop.na.tax <- config$crop.na.tax
mapq.filt   <- config$mapq.filt; if(!all(is.na(mapq.filt))) {
  mapq.filt   <- (as.integer(gsub(':.*', '', config$mapq.filt)) : as.integer(gsub('.*:', '', config$mapq.filt)) )  # default: NA
}
methods.to.use <- unlist(strsplit(config$methods, ';'))
best.mapq      <- as.logical(config$best.mapq)
CIGAR_points   <- separate(data.frame(stringr::str_split_1(config$CIGAR_points, '; ')), col = 1, into = c('key', 'value'), sep = ' = ') %>%
  spread(key, value)
         

#### configs used
message('Configs used: ')
print(t(config))
##


#### Set up future to use multiple cores  with future
plan(multicore, workers = nproc) 


#### RUN MINITAX ON ALIGMNEMENTS ####

#### 1. Best Alignement Method -->> best precision !
if (any(methods.to.use == 'BestAln')) {
  methods  <- 'BestAln'
  taxa.sum.BestAln <- rbindlist(
    future_lapply(paste0(pardir, '/', bamfiles[]),
                  wrap.fun.complete,
                  config=config,  
                  input='bamfile', 
                  #sample=gsub(pattern, '', gsub('.*/', '', bamfiles)),
                  methods=methods, 
                  saveRDS=F, steps=c('1.Import_Alns', '2.Find_taxidentity','3.Refine', '4.Summarise'),
                  outputs=unlist(strsplit(config$outputs, ';')),
                  CIGAR_points=CIGAR_points,
                  best.mapq=T, mapq.filt=mapq.filt,
                  db=db, db.data=db.data, db.uni.data=db.uni.data, 
                  ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    ))
  
  taxa.sum.BestAln$species[is.na(taxa.sum.BestAln$tax.identity)] <- 'unclassified'
  taxa.sum.BestAln$tax.identity[is.na(taxa.sum.BestAln$tax.identity)] <- 'unclassified'
  stopifnot(all(taxa.sum.BestAln$tax.identity == taxa.sum.BestAln$species))
  
  #### Generate phyloseq objects ####
  ps <- PS_from_taxa.sum(taxa.sum.BestAln, metadata = metadata,
                         count.colname  ='count', taxID='species', 
                         cols_to_unique = c(ranks),
                         ranks=ranks)
  sample_data(ps)$method <- methods
  sample_names(ps)       <- paste(sample_data(ps)$workflow, sample_data(ps)$db, sample_data(ps)$method, sample_names(ps), sep='_')
   
  ###
  saveRDS(ps, paste0(outdir, '/', outdir, '_', methods, '_PS.rds'))
  ps.minitax.BA <- ps

}
####


#### 2. LCA Method  -->> ???
if (any(methods.to.use == 'LCA')) {
  ###
  methods  <- 'LCA'
  
  taxa.sum.LCA <- rbindlist(
    future_lapply(paste0(pardir, '/', bamfiles[]),
                  wrap.fun.complete,
                  config=config,  
                  input='bamfile', 
                  #sample=gsub(pattern, '', gsub('.*/', '', bamfiles)),
                  methods=methods, 
                  saveRDS=F, steps=c('1.Import_Alns', '2.Find_taxidentity', '3.Refine', '4.Summarise'),
                  outputs=unlist(strsplit(config$outputs, ';')),
                  CIGAR_points=CIGAR_points,
                  best.mapq=T, mapq.filt=mapq.filt,
                  db=db, db.data=db.data, db.uni.data=db.uni.data, 
                  ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    ))
  taxa.sum.LCA <- rename_duptaxa(taxa.sum.LCA, ranks=ranks)
  
  #### 
  ps <- PS_from_taxa.sum(taxa.sum.LCA, metadata = metadata,
                         count.colname  ='count', taxID='tax.identity', 
                         cols_to_unique = c('tax.identity', ranks),
                         ranks=ranks)
  sample_data(ps)$method <- methods
  sample_names(ps)       <- paste(sample_data(ps)$workflow, sample_data(ps)$db, sample_data(ps)$method, sample_names(ps), sep='_')
  
  saveRDS(ps, paste0(outdir, '/', outdir, '_', methods, '_PS.rds'))
  ps.minitax.LCA <- ps
  
}
####


#### 3. Species Estimation Method -->> much faster !
if (any(methods.to.use == 'SpeciesEstimate')) {

  ###
  methods  <- 'SpeciesEstimate'
  taxa.sum.SpecEst <- rbindlist(
    future_lapply(paste0(pardir, '/', bamfiles[]),
                  wrap.fun.complete,
                  config=config,  
                  input='bamfile', 
                  #sample=gsub(pattern, '', gsub('.*/', '', bamfiles)),
                  methods=methods, 
                  saveRDS=F, steps=c('1.Import_Alns', '2.Find_taxidentity', '3.Refine', '4.Summarise'),
                  outputs=unlist(strsplit(config$outputs, ';')),
                  CIGAR_points=CIGAR_points,
                  best.mapq=T, mapq.filt=mapq.filt,
                  db=db, db.data=db.data, db.uni.data=db.uni.data, 
                  ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
   ))
  taxa.sum.SpecEst$species[is.na(taxa.sum.SpecEst$species)] <- 'unclassified'
  
  
  ps <- PS_from_taxa.sum(taxa.sum.SpecEst, metadata = metadata,
                         count.colname  ='norm_count', taxID='species', 
                         cols_to_unique = c(ranks),
                         ranks=ranks)
  sample_data(ps)$method <- methods
  sample_names(ps)       <- paste(sample_data(ps)$workflow, sample_data(ps)$db, sample_data(ps)$method, sample_names(ps), sep='_')
  
  saveRDS(ps, paste0(outdir, '/', outdir, '_', methods, '_PS.rds'))
  ps.minitax.SE <- ps
  
  
}
####
#### ####


#### Final Output ####
       if (all(c('BestAln', 'LCA', 'SpeciesEstimate')  %in% methods.to.use )) { 
  ps.minitax.comb <- merge_phyloseq(ps.minitax.BA,  ps.minitax.LCA, ps.minitax.SE)
} else if (all(c('BestAln', 'LCA') %in% methods.to.use ))  { 
  ps.minitax.comb <- merge_phyloseq(ps.minitax.BA,  ps.minitax.LCA) 
} else if (all(c('BestAln', 'SpeciesEstimate') %in% methods.to.use )) { 
  ps.minitax.comb <- merge_phyloseq(ps.minitax.BA,  ps.minitax.SE) 
} else if (all(c('LCA', 'SpeciesEstimate') %in% methods.to.use)) { 
  ps.minitax.comb <- merge_phyloseq(ps.minitax.LCA, ps.minitax.SE)
} else if (methods.to.use == 'BestAln') { 
  ps.minitax.comb <- ps.minitax.BA  
} else if (methods.to.use == 'SpeciesEstimate') { 
  ps.minitax.comb <- ps.minitax.SE  
} else if (methods.to.use == 'LCA') {
  ps.minitax.comb <- ps.minitax.LCA 
}

metadata <- data.frame(sample_data(ps.minitax.comb))
write_tsv(metadata, paste0(outdir, "/metadata.tsv"))

if(!all(is.na(mapq.filt)))  {
  save.image(paste0(outdir, '/', outdir, '.MAPQ.filt.RData'))
} else {
  save.image(paste0(outdir, '/', outdir, '.nofilt.RData'))
}

rm(list=grep('ps.minitax.comb|outdir|mapq.filt', ls(), invert = T, value = T))

if(!all(is.na(mapq.filt)))  {
  save.image(paste0(outdir, '/', outdir, '.MAPQ.filt_PS.RData'))
} else {
  save.image(paste0(outdir, '/', outdir, '.nofilt_PS.RData'))
}

################################

