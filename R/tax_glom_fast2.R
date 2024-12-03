tax_glom_fast2 <- function(ps, rank=NULL, rank_level=NULL, ignore_lineage=F) {

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

  ## Lineage, considering all ranks from the top to the supplied rank (level)
  ranks    <- ranks[1:rank_level]

  if (ignore_lineage) {

    # Lineage is the supplied rank
    taxtab[, lineage := get(rank)]
    # Check
    stopifnot(nrow(taxtab[get(rank) != lineage]) == 0)
    # NA, where the lineage is the supplied rank only
    taxtab[is.na(lineage), lineage := 'NA-lineage']

    ranks <- rank

  } else {

    # Create the lineage column by pasting the specified columns using lapply
    #psdata[is.na(get(rank)),  domain := NA]
    taxtab[, lineage := apply(.SD, 1, function(row) paste(na.omit(row), collapse = ";")), .SDcols = ranks]

    # unknown lineage, where the rank above
    taxtab[is.na(get(rank)),  lineage := paste('unknown', lineage, rank, sep = '_')]

  }
  
  
  # merge with counts
  psdata   <- merge(taxtab, otutab, by='rownames', all=T)
  stopifnot(all(!is.na(psdata$rownames)))
  

  # melt for summary
  cols_to_group_by <- c('lineage', ranks)
  psdata.m <- melt.data.table(psdata, id.vars = cols_to_group_by, measure.vars = samples, variable.name = 'sample', value.name = 'count')

  # summarse on given level
  cols_to_sum  <- c(cols_to_group_by, 'sample')
  psdata.m.glom  <- psdata.m[,.(glom_count = sum(count)), by=cols_to_sum]

  # Construct the formula dynamically for Dcast
  formula_str <- paste(paste(cols_to_group_by, collapse = "+"), " ~ sample")
  formula <- as.formula(formula_str)

  psdata.glom <- dcast.data.table(psdata.m.glom, formula, value.var = 'glom_count')

  ## make PS
  taxtab   <- data.frame(psdata.glom[,..ranks  ], row.names = psdata.glom$lineage)
  otutab   <- data.frame(psdata.glom[,..samples], row.names = psdata.glom$lineage)

  colnames(otutab)[]   <- samp_names

  ps.glom <- phyloseq(tax_table(as.matrix(taxtab)),
                      otu_table(as.matrix(otutab), taxa_are_rows = T),
                      sample_data(sampdat))

  return(ps.glom)
}
