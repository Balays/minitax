
### Add taxonomy info (lineage) based on given rank and rank name 

add_lineage <- function(data, ranks, db.uni.data, only.tax.identity=F) {
  
  taxident.all.sum <- data
  
  db.uni  <- unique(db.uni.data[,..ranks]) # c("superkingdom", "phylum", "class", "order", "family", "genus", "species")])
  
  taxident.all.sum.merge <- data.table(NULL)
  for (i in seq_along(ranks)) {
    #rank <- 'species'
    rank <- ranks[i]#'genus'
    message('adding lineage, based on: ', rank, '...')
    db.cols <- ranks[1:i]
    db.uni.rank <- unique(db.uni[,..db.cols])
    
    ## use the rank columns or only the tax.identity?
    if (only.tax.identity) {
      
      cols_by <- 'tax.identity'
      cols    <- colnames(taxident.all.sum)
      cols    <- unique(c(cols, cols_by))
      db.cols <- tail(db.cols, 1)
      
    } else {
    
      cols_by <- c(db.cols[-length(db.cols)], 'tax.identity')
      cols    <- colnames(taxident.all.sum)[!colnames(taxident.all.sum) %in% ranks]
      cols    <- unique(c(cols, cols_by))
      
    }
    
    taxident.sum <- taxident.all.sum[tax.identity.level == rank, ..cols]
    taxident.sum <- merge(db.uni.rank, taxident.sum, by.x=db.cols, by.y=cols_by)
                          #allow.cartesian=TRUE) # 'species'
    taxident.sum[,tax.identity := get(rank)]
    
    if (nrow(taxident.sum) > 0 ) {
      taxident.all.sum.merge <- rbind(taxident.all.sum.merge, taxident.sum, fill=T)
    }
    
  }
  taxident.all.sum.merge[taxident.all.sum.merge == ''] <- NA
  
  taxident.all.sum <- taxident.all.sum.merge
  
  return(taxident.all.sum)
  
}
