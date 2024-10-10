
# Function to add a new column with the rightmost non-NA value
add_rightmost_non_na <- function(dt, columns, new_col_name) {
  dt[, (new_col_name) := apply(.SD, 1, function(row) {
    rev_row <- rev(row)
    rightmost_non_na <- rev_row[which(!is.na(rev_row))[1]]
    return(rightmost_non_na)
  }), .SDcols = columns]
}

# Function to add a new column with the column name of the rightmost non-NA value
add_rightmost_non_na_colname <- function(dt, columns, new_col_name) {
  dt[, (new_col_name) := apply(.SD, 1, function(row) {
    for (i in length(row):1) {
      if (!is.na(row[[i]])) {
        return(columns[i])
      }
    }
    return(NA)
  }), .SDcols = columns]
}
