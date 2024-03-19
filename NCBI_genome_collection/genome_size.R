setwd('E:/data/databases/all_NCBI_genomes')

library(data.table)
library(dplyr)

# Read the contig table
## several 'viruses' were omitted manually !
db.contigs <- fread("NCBI.db.tsv")
ranks <- colnames(db.contigs)[-c(1:4)]
# Read the sequence lengths
## linux: zcat all_NCBI_genomes.fna.gz | bioawk -c fastx '{print $name, length($seq)}' > output.txt
lengths.contigs <- fread('output.txt', col.names = c('seqnames', 'seqlengths'))

db.lengths.contigs <- merge(lengths.contigs, db.contigs, by='seqnames', all=T)

genome_sizes <- #as.data.frame(
  db.lengths.contigs %>% 
  group_by(across(all_of(c(ranks, 'taxid', 'ident')))) %>% 
  summarise(genome_size=sum(seqlengths)) 
#)
## there were some seqnames in the length table which were not present in the DB
genome_sizes <- genome_sizes[!is.na(genome_sizes$ident), ]
genome_sizes[genome_sizes == ''] <- NA

fwrite(genome_sizes, quote = F, sep = '\t', row.names = F,
       'NCBI.db.genomesize.tsv')


