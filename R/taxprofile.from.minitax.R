

# Taxonomic Profiling Output
taxprofile.from.minitax <- function(sample, add.sample.col=F, na.rm=T, 
                                    ranks=c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
                                    SampleID=stri_extract(sample, regex = "\\d+"),
                                    header = data.frame(header=c(
                                      '# Taxonomic Profiling Output',
                                      paste0('@SampleID: ', SampleID),
                                      '@Version:0.9.1',
                                      '@Ranks:superkingdom|phylum|class|order|family|genus|species',
                                      '@TaxonomyID:ncbi-taxonomy_DATE',
                                      '@__program__:minitax',
                                      '@@TAXID RANK TAXPATH TAXPATHSN PERCENTAGE'
                                    )),
                                    outdir    = '.',
                                    infile    = paste0(outdir, "/", sample, ".minitax.tsv"),
                                    outfile   = paste0(outdir, "/", sample, ".tax.profile"),
                                    overwrite = T
                                    ) {
  
  header[7,] <- gsub(' ', '\t', header[7,])
  
  minitax.filt         <- read.delim(infile)
  minitax.uni          <- unique.data.frame(minitax.filt[,ranks])
  
  minitax   <- data.frame(NULL)
  i <- 2
  
  if(add.sample.col) {minitax.filt <- data.frame(sample=sample, minitax.filt)}
  
  for (i in seq_along(ranks)) {
    rank <- ranks[i]
    #minitax.gt         <- minitax.uni  %>% gather(rank, taxon, !all_of(rank))
    minitax.gt         <- unique.data.frame(as.data.frame(minitax.uni[,1:i]))
    colnames(minitax.gt) <- ranks[1:i]
    
    minitax.path       <- as.data.frame(minitax.gt)
    minitax.path       <- as.data.frame(apply(minitax.path, 2, function(x) stri_extract(as.character(unlist(x)), regex = "\\d+")))
    colnames(minitax.path) <- ranks[1:i]
    minitax.path       <- tidyr::unite(minitax.path, 'TAXPATH',       sep = '|', remove = F, na.rm = T, any_of(ranks))
    minitax.path       <- minitax.path[,c(rank, 'TAXPATH')]
    
    minitax.pathsn     <- minitax.gt
    minitax.pathsn     <- as.data.frame(apply(minitax.pathsn, 2, function(x) str_replace(as.character(unlist(x)), pattern = "^\\d+\\s", replacement = "")))
    colnames(minitax.pathsn) <- ranks[1:i]
    minitax.pathsn     <- tidyr::unite(minitax.pathsn, 'TAXPATHSN',   sep = '|', remove = F, na.rm = T, any_of(ranks))
    minitax.pathsn     <- minitax.pathsn[,c(rank, 'TAXPATHSN')]
    
    
    minitax.perc         <- minitax.filt %>% group_by(sample, !!sym(rank)) %>% summarise(count=n(), PERCENTAGE=count/nrow(minitax.filt))
    minitax.perc$TAXID   <- stri_extract(as.character(unlist(minitax.perc[,rank])), regex = "\\d+")
    minitax.perc$RANK    <- rank
    minitax.perc$TAXN    <- str_replace(as.character(unlist(minitax.perc[,rank])), pattern = "^\\d+\\s", replacement = "")
    
    minitax.rank    <- merge(minitax.perc, minitax.path,   by.x='TAXID', by.y=rank)
    minitax.rank    <- merge(minitax.rank,      minitax.pathsn, by.x='TAXN',  by.y=rank)
    minitax.rank    <- minitax.rank[, c('sample', 'TAXN', 'TAXID', 'RANK','TAXPATH', 'TAXPATHSN', 'PERCENTAGE', 'count')]
    
    minitax <- plyr::rbind.fill(minitax, minitax.rank)
  
  }
  
  minitax.taxprofile <- minitax[, c('TAXID', 'RANK','TAXPATH', 'TAXPATHSN', 'PERCENTAGE')]
  
  if (na.rm) {
    minitax.taxprofile <- minitax.taxprofile[!is.na(minitax.taxprofile$TAXID),]
  }
  
  if(overwrite) {
    file.remove(outfile)
  }
  write.table(header, outfile, row.names = F, col.names = F, sep='\t', quote = F)
  write.table(minitax.taxprofile, outfile, append = T, row.names = F, col.names = F, sep='\t', quote = F)
  
}
