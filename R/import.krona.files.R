import.krona.files <- function(
    indir='host_removed_kraken2',
    krona.files = list.files(indir, pattern = '*_krona.txt', full.names = T),
    #lineage.files = list.files(indir, pattern = 'lineage_summary.tsv', full.names = T),
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
  for (i in 1:length(krona.files)) {

    sample <- gsub('.*/', '', gsub('_krona.txt', '', krona.files[i]))

    message('analyzing: ', sample, '...')

    sampdat.sample <- data.frame(sampdat.all,
                                 sample=sample,
                                 row.names=sample)

    data <- fread(krona.files[i], fill=TRUE, na.strings = '')
    colnames(data)[1] <- c('count')

    krona <- unite(as.data.frame(data), 'lineage', -1, sep=';', remove = F, na.rm = T)
    setDT(krona)

    melted <- melt.data.table(krona, id.vars=c('lineage', 'count'), value.name = 'taxon', variable.name = 'rank')
    melted <- melted[!is.na(taxon)]
    melted[, rank := fifelse(grepl('k__', taxon), 'superkingdom',
                             fifelse(grepl('p__', taxon), 'phylum',
                                     fifelse(grepl('c__', taxon), 'class',
                                             fifelse(grepl('o__', taxon), 'order',
                                                     fifelse(grepl('f__', taxon), 'family',
                                                             fifelse(grepl('g__', taxon), 'genus',
                                                                     fifelse(grepl('s__', taxon), 'species',
                                                                             fifelse(grepl('unclassified', taxon, ignore.case = T), 'superkingdom', 'NA'))))))))]

    melted[, taxon := gsub('.__', '', taxon) ]


    spreaded <- dcast.data.table(melted, lineage + count ~ rank, value.var = 'taxon')

    spreaded[grepl('unclassified', superkingdom, ignore.case = T), species := 'Unclassified']

    #spreaded <- spreaded[!is.na(species), ]

    krona <- data.frame(spreaded, row.names = spreaded$lineage)

    ranks <- unique(melted[,rank])

    taxtab <- krona[,ranks]
    otutab <- as.data.frame(krona[,'count']); rownames(otutab) <- rownames(taxtab)
    colnames(otutab) <- sample

    ps <- phyloseq(tax_table(as.matrix(taxtab)),
                   otu_table(as.matrix(otutab),  taxa_are_rows = T),
                   sample_data(sampdat.sample))


    ps.all <- merge_phyloseq(ps, ps.all)

  }

  ps.all <- prune_samples(grep('dummy', sample_names(ps.all), value = T, invert = T), ps.all)

  return(ps.all)
}
