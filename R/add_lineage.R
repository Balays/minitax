
### Add taxonomy info (lineage) based on given rank and rank name

add_lineage <- function(data, ranks, db.uni.data, only.tax.identity = FALSE) {

  taxident.all.sum <- data.table::as.data.table(data.table::copy(data))
  db.uni.data <- data.table::as.data.table(data.table::copy(db.uni.data))

  ## Joins may compare tax.identity with taxid for the taxid rank.
  ## NCBI tables can import taxid as integer while tax.identity is usually character.
  ## Coerce taxonomy/join columns locally to avoid data.table incompatible join-type errors.
  join.candidates <- unique(c("tax.identity", "taxid", ranks))
  for (col in intersect(join.candidates, colnames(taxident.all.sum))) {
    taxident.all.sum[, (col) := as.character(get(col))]
  }
  for (col in intersect(join.candidates, colnames(db.uni.data))) {
    db.uni.data[, (col) := as.character(get(col))]
  }

  db.uni <- unique(db.uni.data[, ..ranks])

  taxident.all.sum.merge <- data.table::data.table(NULL)
  for (i in seq_along(ranks)) {
    rank <- ranks[i]
    message('adding lineage, based on: ', rank, '...')
    db.cols <- ranks[1:i]
    db.uni.rank <- unique(db.uni[, ..db.cols])

    ## use the rank columns or only the tax.identity?
    if (only.tax.identity) {
      cols_by <- 'tax.identity'
      cols <- colnames(taxident.all.sum)
      cols <- unique(c(cols, cols_by))
      db.cols <- tail(db.cols, 1)
    } else {
      cols_by <- c(db.cols[-length(db.cols)], 'tax.identity')
      cols <- colnames(taxident.all.sum)[!colnames(taxident.all.sum) %in% ranks]
      cols <- unique(c(cols, cols_by))
    }

    taxident.sum <- taxident.all.sum[tax.identity.level == rank, ..cols]

    for (col in intersect(db.cols, colnames(db.uni.rank))) {
      db.uni.rank[, (col) := as.character(get(col))]
    }
    for (col in intersect(cols_by, colnames(taxident.sum))) {
      taxident.sum[, (col) := as.character(get(col))]
    }

    taxident.sum <- merge(db.uni.rank, taxident.sum, by.x = db.cols, by.y = cols_by)
    taxident.sum[, tax.identity := get(rank)]

    if (nrow(taxident.sum) > 0) {
      taxident.all.sum.merge <- data.table::rbindlist(
        list(taxident.all.sum.merge, taxident.sum),
        fill = TRUE
      )
    }
  }

  taxident.all.sum.merge[taxident.all.sum.merge == ''] <- NA
  taxident.all.sum.merge
}
