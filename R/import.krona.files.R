import.krona.files <- function(
    indir='host_removed_kraken2',
    krona.files = list.files(indir, pattern = '*_krona.txt', full.names = T),
    #lineage.files = list.files(indir, pattern = 'lineage_summary.tsv', full.names = T),
    ranks = c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
    omit.unclassified = F,
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

    data <- fread(krona.files[i], fill=Inf, na.strings = '')
    colnames(data)[1] <- c('count')

    krona <- unite(as.data.frame(data), 'lineage', -1, sep=';', remove = F, na.rm = T)
    setDT(krona)


    melted <- melt.data.table(krona, id.vars=c('lineage', 'count'), value.name = 'taxon', variable.name = 'rank')
    melted <- melted[!is.na(taxon)]
    melted[, rank_level := as.integer(gsub('V', '', rank))]
    melted[, rank := fifelse(grepl('k__', taxon), 'superkingdom',
                             fifelse(grepl('p__', taxon), 'phylum',
                                     fifelse(grepl('c__', taxon), 'class',
                                             fifelse(grepl('o__', taxon), 'order',
                                                     fifelse(grepl('f__', taxon), 'family',
                                                             fifelse(grepl('g__', taxon), 'genus',
                                                                     fifelse(grepl('s__', taxon), 'species',
                                                                             fifelse(grepl('x__root', taxon, ignore.case = T), 'root',
                                                                                     fifelse(grepl('x__', taxon, ignore.case = T), 'other',
                                                                                             fifelse(grepl('unclassified', taxon, ignore.case = T), ranks[1],
                                                                                     'NA'))))))))))]

    melted[, species_rank := max(fifelse(rank == "species", rank_level, NA_integer_), na.rm = TRUE), by = lineage]
    melted[is.infinite(species_rank), species_rank := NA_integer_]  # if a lineage has no species
    melted[, has_species := !is.na(species_rank)]
    melted[, has_subspecies := has_species & any(rank_level == species_rank + 1L & rank == "other"), by = lineage]
    melted[, subspecies_rank := fifelse(has_subspecies, species_rank + 1L, NA_integer_), by = lineage]

    ## make subsepcies those that have a rank level higher than species
    melted[
      !is.na(subspecies_rank) &
        rank == "other" &
        rank_level == subspecies_rank,
      rank := "subspecies"
    ]


    ###
    #melted[, rank_levels      := paste(rank, collapse = '::'), by = lineage ]
    #melted[, has_species      := fifelse(grepl('species', rank_levels), T, F), by = lineage]
    #melted[, has_subspecies   := fifelse(grepl('species::other', rank_levels), T, F), by = lineage]
    #melted[, species_rank     := fifelse(rank == 'species', rank_level,  0), by = lineage]
    #melted[, subspecies_rank  := fifelse(has_subspecies == T, species_rank  + 1, 0), by = lineage]
    #taxa_dt <- unique(melted[,.(rank, taxon)])

    ##
    melted[, taxon := gsub('.__', '', taxon) ]

    ## Keep only 'ranks'
    melted   <- melted[rank %in% ranks, ]

    melted[, rank := factor(rank, levels = ranks)]

    melted   <- melted[order(lineage, rank)]

    #melted[, lineage := paste0( sep=';')]

    spreaded <- dcast.data.table(melted, lineage + count ~ rank, value.var = 'taxon')

    spreaded[grepl('unclassified', ranks[1], ignore.case = T), species := 'Unclassified']

    # Filtering
    #spreaded <- spreaded[!is.na(species), ]

    if(exists('omit.unclassified') && omit.unclassified) {
      spreaded <- spreaded[!grepl('unclassified', lineage, ignore.case = T), ]
    }


    melted[, rank := as.character(rank)]
    ranks <- unique(melted[,rank])
    cols  <- c('lineage', 'count', as.character(ranks))

    krona <- unique(spreaded[,..cols])

    krona <- setDT(unite(data.frame(krona), col = 'lineage', all_of(ranks), sep=';', remove = F, na.rm = T))

    krona <- unique(krona[,..cols])

    if(length(dup(krona$lineage)) > 0) {
      message('There are ', length(dup(krona$lineage)), ' taxa, where there is a difference in the lineages, other than specified in the "ranks" argument.\n
              These are needed to be collapsed, so it was carried out.')

      cols  <- c('lineage', as.character(ranks))

      counts <- krona[, .(count = sum(count)), by=cols ]

      krona  <- counts
    }

    krona <- data.frame(krona, row.names = krona$lineage)


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
