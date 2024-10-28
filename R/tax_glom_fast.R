

tax_glom_fast <- function(ps, rank=NULL, rank_level=NULL, ignore_lineage=T) {
  
  require(data.table)
  
  ranks    <- rank_names(ps)
  
  if (is.null(rank) && is.null(rank_level)) {
    stop("Either 'rank' or 'rank_level' must be provided.")
  }
  
  if (!is.null(rank) && !is.null(rank_level)) {
    message("Both 'rank' and 'rank_level' are provided. Using 'rank' to determine 'rank_level'.")
  }
  
  if (!is.null(rank)) {
    rank_level <- which(ranks == rank)
    if (length(rank_level) == 0) {
      stop("The provided 'rank' does not exist in 'ranks'.")
    }
  } else {
    if (rank_level > length(ranks) || rank_level < 1) {
      stop("The provided 'rank_level' is out of range.")
    }
    rank <- ranks[rank_level]
  }
  
  samp_names         <- sample_names(ps)
  
  otutab   <- data.frame(otu_table(ps)); otutab   <- data.table(rownames=rownames(otutab), otutab)
  taxtab   <- data.frame(tax_table(ps)); taxtab   <- data.table(rownames=rownames(taxtab), taxtab)
  sampdat  <- data.frame(sample_data(ps))
  samples  <- rownames(sampdat)
  
  colnames(otutab)[-1]   <- samp_names
  
  psdata   <- merge(taxtab, otutab, by='rownames', all=T)
  stopifnot(all(!is.na(psdata$rownames)))
  
  ## Lineage, considering all ranks from the top to the supplied rank (level)
  ranks    <- ranks[1:rank_level]
  
  if (!ignore_lineage) {
    # Create the lineage column by pasting the specified columns using lapply
    psdata[, lineage := apply(.SD, 1, function(row) paste(na.omit(row), collapse = ";")), .SDcols = ranks]
  } else {
    #psdata[, lineage := get(rank)]
    psdata[, lineage := apply(.SD, 1, function(row) paste(na.omit(row), collapse = ";")), .SDcols = ranks]
  }  
  psdata[lineage == '' & is.na(get(rank)),  lineage := 'NA-lineage']
  #psdata[is.na(lineage), lineage := 'NA-lineage']
  
  # melt for summary
  cols_to_group_by <- c('lineage', ranks)
  psdata.m <- melt(psdata, id.vars = cols_to_group_by, measure.vars = samples, variable.name = 'sample', value.name = 'count')
  
  # summarse on given level
  cols_to_sum  <- c(cols_to_group_by, 'sample')
  psdata.m.glom  <- psdata.m[,.(glom_count = sum(count)), by=cols_to_sum]
  
  # Construct the formula dynamically for Dcast
  formula_str <- paste(paste(cols_to_group_by, collapse = "+"), " ~ sample")
  formula <- as.formula(formula_str)
  
  psdata.glom <- dcast(psdata.m.glom, formula, value.var = 'glom_count')
  
  ## make PS
  taxtab   <- data.frame(psdata.glom[,..ranks  ], row.names = psdata.glom$lineage)
  otutab   <- data.frame(psdata.glom[,..samples], row.names = psdata.glom$lineage)
  
  colnames(otutab)[]   <- samp_names
  
  ps.glom <- phyloseq(tax_table(as.matrix(taxtab)),
                      otu_table(as.matrix(otutab), taxa_are_rows = T),
                      sample_data(sampdat))
  
  return(ps.glom)
}