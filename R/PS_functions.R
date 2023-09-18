
## merged dataframe fro phyloseq object
df.from.ps <- function(ps,
                       #metadata=NA,
                       by='sample',
                       n=NA, all=T,
                       comb=T, i=NA, t=NA, NArm = FALSE) {

  if(is.na(n))      { n      <- length(rank_names(ps)) }

  if(comb) {
    if(is.na(i))   { i  <- 2      }
    if(is.na(t))   { t  <- rank_names(ps)[i]      }
    ps <- tax_glom(ps, t, NArm = NArm)
    message('Species were combined at ', t, ' level.')
  }

  taxtable <- as.data.frame(tax_table(ps)[,c(1:n)])
  otutable <- as.data.frame(otu_table(ps, taxa_are_rows = T))

  psdf <- merge(
    taxtable,
    otutable,
    by=0, all=all
  )[,-1]

  return(psdf)
}



## normalize counts
norm.ps <- function(ps, norm=NA, pseudocount=0) {
  otutab <- otu_table(ps) + pseudocount
  ps <- phyloseq(otutab, tax_table(ps), sample_data(ps))
  if (is.na(norm)) {
    ps.prop  <- ps
  } else if (norm == 'sum') {
    ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
  } else if (norm == 'rlog') {
    ps.prop  <- ps
    de.seq   <- phyloseq_to_deseq2(ps, ~sample)
    rlog.asv <- assay(rlog(de.seq))
    otu_table(ps.prop) <- otu_table(rlog.asv, taxa_are_rows = T)
  } else  {stop('incorrect "norm"!')}

  return(ps.prop)
}


colordf <- data.frame(npg        = pal_npg()(10),
                      nejm       = pal_nejm()(10),
                      aas        = pal_aaas()(10),
                      jco        = pal_jco()(10),
                      uchicago   = pal_uchicago()(10),
                      locuszoom  = pal_locuszoom()(10),
                      lancet     = pal_lancet()(10)
                      )

colorvec <- as.character(unlist(lapply(colordf, c)))
#scales::show_col(colorvec)



PS_from_taxa.sum <- function(
    taxa.sum, metadata=NULL,
    count.colname  ='count', ranks=NULL,
    cols_to_unique = c('tax.identity', ranks),
    taxID='species',
    round_counts_decimal=0
) {

  require(phyloseq)

  ### Taxonomy table ####
  taxtab   <- unique(taxa.sum[,..cols_to_unique])
  taxtab   <- data.frame(taxtab[,], row.names = unlist(taxtab[,..taxID]))


  ### Count (OTU) table ####

  if (!is.na(round_counts_decimal)) {
    taxa.sum[, (count.colname) := lapply(.SD, round, round_counts_decimal), .SDcols = count.colname]
  }

  cols_to_sum <- c(taxID, 'sample', count.colname)
  otutab    <- data.table(taxa.sum[,..cols_to_sum] %>% spread('sample',count.colname, fill=0))
  otutab    <- data.frame(otutab[,-1], row.names = unlist(otutab[,..taxID]))

  if(length(samples) == 1 ) {colnames(otutab) <- samples}

  ### Metadata ####
  if(is.null(metadata)) {
    metadata <- data.frame(sample=unique(taxa.sum$sample), row.names = unique(taxa.sum$sample))
  }

  ### make PS object ####

  ps <- phyloseq(tax_table(as.matrix(taxtab)),
                 otu_table(as.matrix(otutab), taxa_are_rows = T),
                 sample_data(metadata))


  return(ps)
}


## function to glomerate a PS object at n taxonomic level, optionally regardless of lineage
glom_taxa_n <- function(ps, n=2, ignore_lineage=T) {

  t <- rank_names(ps)[n]
  ps.glom   <- tax_glom(ps, t, NArm = F)

  taxtab.glom    <- as.data.frame(tax_table(ps.glom)[,c(1:n)])
  otutab.glom    <- as.data.frame(otu_table(ps.glom)[,])

  ## change rownames to a unique lineage (which HAS to be unique since we used tax_glom )
  taxtab.glom <- unite(taxtab.glom, 'lineage', sep='::', remove = F, na.rm = T)
  # use existing rownames if all taxon levels are NA
  taxtab.glom$lineage[taxtab.glom$lineage == ''] <- rownames(taxtab.glom)[taxtab.glom$lineage == '']
  # keep the order, just in case
  otutab.glom <- as.data.frame(otutab.glom[rownames(taxtab.glom),])
  rownames(otutab.glom) <- rownames(taxtab.glom)
  colnames(otutab.glom) <- sample_names(ps.glom)

  ## rename to the specific tax level
  rownames(taxtab.glom) <- taxtab.glom[,'lineage']
  rownames(otutab.glom) <- taxtab.glom[,'lineage']
  taxtab.glom <- dplyr::select(taxtab.glom, -lineage)

  ## regenerate PS object
  ps.glom <- phyloseq(tax_table(as.matrix(taxtab.glom)),
                      otu_table(as.matrix(otutab.glom), taxa_are_rows = T),
                      sample_data(ps.glom))

  ## tax glom again, without higher taxon levels in order to glom together the same taxa with differing lineages
  if(ignore_lineage) {
    taxtab.glom       <- data.frame(dummy='', tax_table(ps.glom)[,n])
    ps.glom@tax_table <- tax_table(as.matrix(taxtab.glom))

    taxtab.glom       <- as.data.frame(tax_table(ps.glom)[,])
    otutab.glom       <- as.data.frame(otu_table(ps.glom)[,])

    ps.glom.glom      <- tax_glom(ps.glom, t, NArm = F)
    taxtab.glom.glom  <- as.data.frame(tax_table(ps.glom.glom)[,])
    otutab.glom.glom  <- as.data.frame(otu_table(ps.glom.glom)[,])

    ## change rownames to a unique lineage (which HAS to be unique since we used tax_glom() )
    taxtab.glom.glom <- unite(taxtab.glom.glom %>% dplyr::select(-dummy), 'lineage', sep='::', remove = F, na.rm = T)

    ## keep the order, just in case
    otutab.glom.glom <- as.data.frame(otutab.glom.glom[rownames(taxtab.glom.glom),])
    rownames(otutab.glom.glom) <- rownames(taxtab.glom.glom)
    colnames(otutab.glom.glom) <- sample_names(ps.glom.glom)


    ## rename to the specific tax level
    rownames(taxtab.glom.glom) <- taxtab.glom.glom[,'lineage']
    rownames(otutab.glom.glom) <- taxtab.glom.glom[,'lineage']
    taxtab.glom.glom <- dplyr::select(taxtab.glom.glom, -lineage)
    ## re-regenerate PS object
    ps.glom.glom <- phyloseq(tax_table(as.matrix(taxtab.glom.glom)),
                             otu_table(as.matrix(otutab.glom.glom), taxa_are_rows = T),
                             sample_data(ps.glom.glom))

    ## OKAY !
    ps.glom <- ps.glom.glom
  }

  return(ps.glom)
}



### function to combine samples in a PS object base on one or more variables (columns in sample_data)

glom_samples_var <- function(ps, cols_to_group_by, by.x='sample', by.y='sample', sum.fun='sum') {


  taxtab    <- data.frame(tax_table(ps))
  otutab    <- data.frame(taxon=taxa_names(ps), otu_table(ps))
  colnames(otutab)[-1] <- sample_names(ps)
  sampdat   <- data.frame(sample_data(ps))

  otutab.gt <- otutab %>% gather(sample, count, -1)
  otutab.gt <- merge(otutab.gt, sampdat, by.x=by.x, by.y=by.y)

  otutab.sum <- unite(otutab.gt, col = 'sample', cols_to_group_by, remove = F, na.rm = T, sep = '_')

  if (sum.fun == 'sum') {
    otutab.sum <- otutab.sum %>% group_by(all_of(across(c('taxon', 'sample', cols_to_group_by)))) %>% summarise(count=sum(count))
  } else if (sum.fun == 'mean') {
    otutab.sum <- otutab.sum %>% group_by(all_of(across(c('taxon', 'sample', cols_to_group_by)))) %>% summarise(count=mean(count))
  }

  sampdat.sum <- as.data.frame(unique.data.frame(otutab.sum[, c('sample', cols_to_group_by)]))
  rownames(sampdat.sum) <- sampdat.sum$sample

  otutab.sum <- unique.data.frame(otutab.sum[, c('taxon', 'sample', 'count')])
  otutab.sum <- spread(otutab.sum, sample, count)
  otutab.sum <- data.frame(otutab.sum[,-1], row.names = otutab.sum$taxon)

  ps.sum <- phyloseq(otu_table(as.matrix(otutab.sum), taxa_are_rows = T),
                     tax_table(as.matrix(taxtab)),
                     sample_data(sampdat.sum)
                     )
  return(ps.sum)
}


### function that glomerates a set of taxa only in PS object at specific rank

glom_sub_taxa <- function(ps, taxa_to_glom, n, ignore_lineage=T) {

  t <- rank_names(ps)[n]
  taxtab           <- data.frame(tax_table(ps))
  species_to_glom  <- rownames(taxtab[is.element(taxtab[,t], taxa_to_glom), ])

  ps.toglom        <- prune_taxa(species_to_glom, ps)
  ps.glom          <- glom_taxa_n(ps.toglom, n, ignore_lineage)
  taxtab           <- data.frame(tax_table(ps.glom))

  taxa_not_glom    <- taxa_names(ps)[!is.element(taxa_names(ps), species_to_glom)]
  ps.ori           <- prune_taxa(taxa_not_glom, ps)

  ps <- merge_phyloseq(ps.ori, ps.glom)

}













































