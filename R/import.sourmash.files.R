import.sourmash.files <- function(
    indir='sourmash_NCBI',
    kreport.files = list.files(indir, pattern = '*.kreport.txt', full.names = T),
    lineage.files = list.files(indir, pattern = 'lineage_summary.tsv', full.names = T),
    ranks = c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
    sampdat.all = data.frame(platform='ONT',
                             Vregion='WGS',
                             workflow='sourmash',
                             method='k31',
                             db='all_NCBI_genomes')

) {

  #read.delim('sourmash_NCBI_k31/ERR5396170.gather.k31.lineage_summary.tsv', header = T)

  ps.all <- phyloseq(tax_table(as.matrix(data.frame(superkingdom = 'Bacteria'))),
                     otu_table(as.matrix(data.frame(dummy = 0)), taxa_are_rows = T))

  i <- 1
  for (i in 1:length(lineage.files)) {

    sample <- gsub('.*/', '', gsub('.kreport.txt', '', kreport.files[i]))

    message('analyzing: ', sample, '...')

    sampdat.sample <- data.frame(sampdat.all,
                                 sample=sample,
                                 row.names=sample)

    #readcounts <- read.delim('readcount_portik.txt')


    data <- fread(kreport.files[i], na.strings = '', header = F)
    colnames(data) <- c('ratio', 'sum.count', 'count', 'rank', 'taxid', 'taxon')

    kreport <- data %>% select(taxon, sum.count, rank) %>% spread(rank, sum.count)
    kreport$sum.count <- rowSums(kreport[,-1], na.rm = T)
    # apply(ps.data[,-1], 1, function(x) sum(x[!is.na(x)], na.rm=T))

    #data$ratio <- round(data$ratio * readcounts$read_count / 100, 0)
    #colnames(kreport)[length(kreport)] <- 'count'

    sum(kreport$count[!is.na(kreport$S)])
    #readcounts$read_count


    ## OK!

    colnames(kreport)[length(kreport)] <- 'kreport'


    sample <- gsub('.*/', '', gsub('.lineage_summary.tsv', '', lineage.files[i]))
    ###

    data <- fread(lineage.files[i], header = T, na.strings = '')
    colnames(data) <- c('lineage', 'fraction')

    lineage_sum <- separate(data, 'lineage', into = ranks, sep=';', remove = T)
    lineage_sum[lineage_sum == ''] <- NA
    #data$count <- round(data$fraction * readcounts$read_count, 0)

    #sum(data$count)
    #readcounts$read_count
    ## OK!

    colnames(lineage_sum)[length(lineage_sum)] <- 'lineage_summary'


    lineage_sum[lineage_sum$superkingdom == 'unclassified', ranks] <- 'unclassified'

    merged.data <- merge(kreport[,c('taxon', 'kreport')], lineage_sum[,], by.x='taxon', by.y='species')
    rownames(merged.data) <- merged.data$taxon
    colnames(merged.data)[1] <- 'species'
    ###

    taxtab.pb.sm <- merged.data[,ranks]
    otutab.pb.sm <- data.frame(count=merged.data[,'kreport'], row.names = rownames(merged.data))
    colnames(otutab.pb.sm) <- sample

    ps.portik.ncbi <- phyloseq(tax_table(as.matrix(taxtab.pb.sm)),
                               otu_table(as.matrix(otutab.pb.sm), taxa_are_rows = T),
                               sample_data(sampdat.sample))


    ps.all <- merge_phyloseq(ps.portik.ncbi, ps.all)

  }

  ps.all <- prune_samples(grep('dummy', sample_names(ps.all), value = T, invert = T), ps.all)

  return(ps.all)
}
