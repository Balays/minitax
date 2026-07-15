


# Function to calculate and assign the most probable taxon for a given rank
assign_most_probable <- function(data, rank, thresholds, ranks) {

  # Define the taxonomic ranks in order
  rank_df <- data.frame(rank=ranks, rank_level = 1:length(ranks))

  # Calculate the total number of alignments for each qname
  total_alignments <- data[, .N, by = .(sample, qname, tax.identity, tax.identity.level, tax.identity_rank_level) ]
  setnames(total_alignments, 'N', 'N_alignments')

  # Calculate the number of alignments for each taxon at the given rank within each qname
  cols_by <- c('sample', 'qname', rank)
  taxon_alignments <- data[!is.na(get(rank)) & get(rank) != '', .(taxon_count = .N), by = cols_by]
  if (nrow(taxon_alignments) == 0) {
    empty <- data.table(
      sample = character(),
      qname = character(),
      tax.identity = character(),
      tax.identity.level = character(),
      tax.identity_rank_level = integer()
    )
    empty[, (rank) := character()]
    empty[, (paste0(rank, '_count')) := integer()]
    empty[, (paste0(rank, '_prob')) := numeric()]
    empty[, (paste0(rank, '_most_probable')) := character()]
    empty[, (paste0(rank, '_most_probable_rank_level')) := integer()]
    empty[, (paste0('is_', rank, '_most_prob')) := logical()]
    return(empty)
  }

  # Merge the total alignments with taxon alignments
  taxon_alignments <- merge(taxon_alignments, total_alignments, by = c("sample", "qname"))

  # Calculate the probability for each taxon at the given rank
  taxon_alignments[, taxon_prob := taxon_count / N_alignments]

  # Assign the most probable taxon to each qname
  taxon_alignments[, taxon_most_probable := get(rank)[which.max(taxon_prob)], by = .(sample, qname)]

  # Assign rank level
  taxon_alignments[,taxon_most_probable_rank := rank]

  taxon_alignments <- merge(taxon_alignments, rank_df, by.x='taxon_most_probable_rank', by.y='rank')
  setnames(taxon_alignments, 'rank_level', "taxon_most_probable_rank_level")

  # Check if the most probable taxon is a better hit than the tax identity
  rank_threshold <- thresholds$threshold[thresholds$rank == rank]
  if (length(rank_threshold) != 1 || is.na(rank_threshold)) {
    stop('Missing BestAln threshold for rank: ', rank)
  }
  taxon_alignments[,is_taxon_most_prob := fifelse(tax.identity_rank_level < taxon_most_probable_rank_level &
                                                    get(rank) == taxon_most_probable &
                                                    taxon_prob > rank_threshold,
                                                  T, F)]

  # Rename the columns appropriately
  setnames(
    taxon_alignments,
    old = c('taxon_count', 'taxon_prob', 'taxon_most_probable',
            'taxon_most_probable_rank', 'taxon_most_probable_rank_level',
            'is_taxon_most_prob'),
    new = c(paste0(rank, '_count'), paste0(rank, '_prob'), paste0(rank, '_most_probable'),
            paste0(rank, '_most_probable_rank'), paste0(rank, '_most_probable_rank_level'),
            paste0('is_', rank, '_most_prob'))
  )

  cols_by <- c('sample', 'qname', 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level',
               grep(rank, colnames(taxon_alignments), value = T))
  return(taxon_alignments[, ..cols_by])
}


# Function to update the tax.identity based on probability and level
update_tax_identity <- function(data, rank) {

  prob_col <- paste0(rank, "_prob")
  most_probable_col <- paste0(rank, "_most_probable")
  most_probable_rank_level_col <- paste0(rank, "_most_probable_rank_level")
  is_taxon_most_prob_col <- paste0('is_', rank, "_most_prob")

  cols <- c('sample', 'qname', 'aln_nr', 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level', rank,
            prob_col, most_probable_col, most_probable_rank_level_col, is_taxon_most_prob_col)

  data.subs <- data[,..cols]

  setnames(data.subs, names(data.subs)[-c(1:6)], gsub(rank, 'taxa', names(data.subs)[-c(1:6)]))

  data.subs[, update_tax.ident := any(is_taxa_most_prob, na.rm = TRUE), by=.(sample, qname)]
  data.subs[, chosen_taxa := {
    hits <- taxa_most_probable[is_taxa_most_prob == TRUE & !is.na(is_taxa_most_prob)]
    if (length(hits) > 0) hits[1] else NA_character_
  }, by=.(sample, qname)]
  data.subs[, chosen_taxa_rank_level := {
    hits <- taxa_most_probable_rank_level[is_taxa_most_prob == TRUE & !is.na(is_taxa_most_prob)]
    if (length(hits) > 0) hits[1] else NA_real_
  }, by=.(sample, qname)]

  data.subs[, tax.identity            := fifelse(update_tax.ident, chosen_taxa, tax.identity)]
  data.subs[, tax.identity.level      := fifelse(update_tax.ident, rank, tax.identity.level)]
  data.subs[, tax.identity_rank_level := fifelse(update_tax.ident, chosen_taxa_rank_level, tax.identity_rank_level)]

  data[,tax.identity             := NULL]
  data[,tax.identity.level       := NULL]
  data[,tax.identity_rank_level  := NULL]

  cols <- c('sample', 'qname', 'aln_nr', 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level')
  data.subs <- data.subs[,..cols]

  data <- merge(data[,], data.subs[,],
                by = c('sample', 'qname', 'aln_nr'))

  return(data)
}


BestAln <- function(taxa.dup, ranks=ranks, thresholds=thresholds) {

  # Define the taxonomic ranks in order
  rank_df <- data.frame(rank=ranks, rank_level = 1:length(ranks))

  taxa.dup     <- merge(taxa.dup, rank_df, by.x='tax.identity.level', by.y='rank', all.x=T)
  setnames(taxa.dup, 'rank_level', 'tax.identity_rank_level')

  taxa.dup[is.na(tax.identity_rank_level), tax.identity_rank_level := 0]

  cols_by <- c('sample', 'qname', 'aln_nr', ranks, 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level')
  taxa.most_prob <- taxa.dup[, ..cols_by]

  ## 1.
  # Apply a Function that calculate and assign the most probable taxon for a given rank, iteratively for each rank
  for (rank in ranks) {
    message('Finding probability for ', rank, '...')
    most_probable_taxa <- assign_most_probable(data=taxa.most_prob, rank, thresholds, ranks)
    taxa.most_prob     <- merge(taxa.most_prob, most_probable_taxa,
                                by = c("sample", "qname", 'tax.identity', 'tax.identity.level', 'tax.identity_rank_level', rank), all.x = TRUE)

  }

  ## 2.
  # Apply a Function to update the tax.identity based on probability and level iteratively for each level
  taxa.new.taxident <- data.table(as.data.frame(taxa.most_prob))

  for (rank in ranks) {
    message('Refining tax.identity based on ', rank, '...')

    taxa.new.taxident <- update_tax_identity(data=taxa.new.taxident, rank)

  }

  return(taxa.new.taxident)

}

minitax_add_taxid_to_summary <- function(taxa.sum, db.uni.data = NULL, ranks = NULL) {
  taxa.sum <- data.table::as.data.table(data.table::copy(taxa.sum))
  if (nrow(taxa.sum) == 0L) return(taxa.sum)

  if (!'taxid' %in% colnames(taxa.sum)) {
    taxa.sum[, taxid := NA_character_]
  }
  taxa.sum[, taxid := as.character(taxid)]

  if (is.null(db.uni.data) || is.null(ranks) || !'taxid' %in% colnames(db.uni.data)) {
    return(taxa.sum)
  }

  db.uni.data <- data.table::as.data.table(data.table::copy(db.uni.data))
  ranks <- intersect(ranks, colnames(db.uni.data))
  if (length(ranks) == 0L) return(taxa.sum)

  ## SpeciesEstimate summaries initially contain only rank columns. Add temporary
  ## tax.identity columns so the same unambiguous-rank mapping can be used there.
  temporary.tax.identity <- FALSE
  if (!all(c('tax.identity', 'tax.identity.level') %in% colnames(taxa.sum))) {
    terminal.rank <- tail(intersect(ranks, colnames(taxa.sum)), 1)
    if (length(terminal.rank) == 1L) {
      taxa.sum[, tax.identity := as.character(get(terminal.rank))]
      taxa.sum[, tax.identity.level := terminal.rank]
      temporary.tax.identity <- TRUE
    }
  }

  if (!all(c('tax.identity', 'tax.identity.level') %in% colnames(taxa.sum))) {
    return(taxa.sum)
  }

  db.uni.data[, taxid := as.character(taxid)]
  rank.maps <- lapply(ranks, function(rank) {
    if (!rank %in% colnames(db.uni.data)) return(data.table::data.table())
    x <- db.uni.data[!is.na(get(rank)) & get(rank) != '' &
                       !is.na(taxid) & taxid != '',
                     .(tax.identity = as.character(get(rank)), taxid = as.character(taxid))]
    x <- unique(x)
    if (nrow(x) == 0L) return(data.table::data.table())
    x <- x[, .(taxid.mapped = if (data.table::uniqueN(taxid) == 1L) unique(taxid) else NA_character_),
           by = tax.identity]
    x[, tax.identity.level := rank]
    x
  })
  taxid.map <- data.table::rbindlist(rank.maps, fill = TRUE)
  taxid.map <- taxid.map[!is.na(taxid.mapped) & taxid.mapped != '']
  if (nrow(taxid.map) == 0L) {
    if (temporary.tax.identity) {
      taxa.sum[, tax.identity := NULL]
      taxa.sum[, tax.identity.level := NULL]
    }
    return(taxa.sum)
  }

  taxa.sum[, row_id__minitax_taxid := .I]
  taxa.sum[, tax.identity := as.character(tax.identity)]
  taxa.sum[, tax.identity.level := as.character(tax.identity.level)]
  taxid.map[, tax.identity := as.character(tax.identity)]
  taxid.map[, tax.identity.level := as.character(tax.identity.level)]

  taxa.sum <- merge(
    taxa.sum,
    taxid.map,
    by = c('tax.identity.level', 'tax.identity'),
    all.x = TRUE,
    sort = FALSE
  )
  taxa.sum[is.na(taxid) | taxid == '', taxid := taxid.mapped]
  taxa.sum[, taxid.mapped := NULL]
  data.table::setorder(taxa.sum, row_id__minitax_taxid)
  taxa.sum[, row_id__minitax_taxid := NULL]

  if (temporary.tax.identity) {
    taxa.sum[, tax.identity := NULL]
    taxa.sum[, tax.identity.level := NULL]
  }

  preferred <- intersect(c('sample', 'lineage', 'taxid'), colnames(taxa.sum))
  data.table::setcolorder(
    taxa.sum,
    unique(c(preferred, setdiff(colnames(taxa.sum), preferred)))
  )
  taxa.sum
}

## Keep large multi-sample runs fault-tolerant during cached summarisation.
## minitax.complete.R sources this file after R/minitax.wrapfun.complete.R, so this
## wrapper catches one failed cached sample without aborting the full future_lapply()
## method run. The failed sample is recorded in classification_timing.tsv with
## status='failed' and the error message; successful samples are still combined.
if (exists('summarise_minitax_taxa_file', mode = 'function') &&
    !exists('.minitax_summarise_minitax_taxa_file_unsafe', inherits = FALSE)) {
  .minitax_summarise_minitax_taxa_file_unsafe <- summarise_minitax_taxa_file

  summarise_minitax_taxa_file <- function(taxa_or_bamfile, ...) {
    args <- list(...)
    methods <- if ('methods' %in% names(args)) args$methods else NA_character_
    method.label <- if (is.null(methods) || length(methods) == 0 || all(is.na(methods))) {
      NA_character_
    } else {
      paste(as.character(methods), collapse = ';')
    }
    sample <- gsub('_best_alignments_w_taxa.tsv$', '', basename(taxa_or_bamfile))
    started <- Sys.time()
    taxa.input <- taxa_or_bamfile
    fixed.cache.dir <- NULL

    tryCatch(
      {
        ## Older cache files may lack a sample column because wrap.fun.complete()
        ## historically left taxa$sample <- sample commented out for input='taxa.DT'.
        ## If such a cache is summarised directly, data.table may pick up base::sample
        ## as a closure in grouping expressions, which later breaks fwrite().
        cache.header <- names(data.table::fread(taxa_or_bamfile, nrows = 0, showProgress = FALSE))
        if (!'sample' %in% cache.header) {
          message('Cached taxa file has no sample column; adding sample from filename for summarisation: ', sample)
          fixed.cache.dir <- tempfile('minitax_cache_sample_')
          dir.create(fixed.cache.dir, recursive = TRUE, showWarnings = FALSE)
          taxa.input <- file.path(fixed.cache.dir, basename(taxa_or_bamfile))
          taxa.cache <- data.table::fread(taxa_or_bamfile, na.strings = '', showProgress = FALSE)
          taxa.cache[taxa.cache == ''] <- NA
          taxa.cache[, sample := sample]
          data.table::fwrite(taxa.cache, taxa.input, sep = '\t', na = 'NA')
          rm(taxa.cache)
          on.exit(unlink(fixed.cache.dir, recursive = TRUE, force = TRUE), add = TRUE)
        }

        result <- .minitax_summarise_minitax_taxa_file_unsafe(taxa.input, ...)
        if (!is.null(result$taxa.sum) && nrow(result$taxa.sum) > 0L) {
          db.uni.data <- if ('db.uni.data' %in% names(args)) args$db.uni.data else NULL
          ranks <- if ('ranks' %in% names(args)) args$ranks else NULL
          result$taxa.sum <- minitax_add_taxid_to_summary(result$taxa.sum, db.uni.data, ranks)
        }
        result
      },
      error = function(e) {
        finished <- Sys.time()
        msg <- conditionMessage(e)
        warning('minitax summarisation failed for ', sample,
                ' [', method.label, ']: ', msg, call. = FALSE)
        list(
          taxa.sum = data.table(),
          timing = data.table(
            sample = sample,
            stage = 'summarise',
            method = method.label,
            status = 'failed',
            error = msg,
            elapsed_sec = as.numeric(difftime(finished, started, units = 'secs')),
            rows_in = NA_integer_,
            rows_out = 0L,
            file = taxa_or_bamfile
          )
        )
      }
    )
  }
}
