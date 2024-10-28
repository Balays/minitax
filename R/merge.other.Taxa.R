

### Merge every Taxa other than those supplied in taxa.to.keep into one

merge.other.Taxa <- function(ps, taxa.to.keep=NULL, other.name='other') {
  
  if (!is.null(taxa.to.keep)) {
    
    try({
      ps.filt <- prune_taxa(taxa.to.keep, ps)
    
      taxa.to.glom <- taxa_names(ps)[!taxa_names(ps) %in% taxa.to.keep]
      ps.other     <- prune_taxa(taxa.to.glom, ps)
      
      taxtab.other          <- data.frame(ps.other@tax_table)
      taxtab.other[,'glom'] <- other.name
      ps.other@tax_table    <- tax_table(as.matrix(taxtab.other))
      
      ps.other.glom   <- tax_glom_fast(ps.other, rank='glom', ignore_lineage = T)
      
      colnames(ps.other.glom@tax_table) <- colnames(ps.filt@tax_table)
      #ps.other.glom@tax_table
      
      ps <- merge_phyloseq(ps.filt, ps.other.glom)
      #rowSums(ps@otu_table)
    })
    
  } else {
    NULL
  }
  return(ps)
}

