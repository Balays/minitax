


# Function to paste together columns, omitting NAs, and insert at a specified position
require(data.table)

# Function to paste together columns, omitting NAs, and insert at a specified position
paste_columns_dt <- function(dt, columns, separator = ",", new_col_name = "pasted_col", position = ncol(dt) + 1) {
  # Create the new column
  dt[, (new_col_name) := apply(.SD, 1, function(row) {
    paste(na.omit(row), collapse = separator)
  }), .SDcols = columns]
  
  # Calculate the correct position
  if (position > ncol(dt)) {
    position <- ncol(dt)
  }
  
  # Reorder columns to place the new column at the specified position
  if (position == 1) {
    setcolorder(dt, c(new_col_name, names(dt)[-ncol(dt)]))
  } else {
    setcolorder(dt, c(names(dt)[1:(position-1)], new_col_name, names(dt)[position:(ncol(dt)-1)]))
  }
}

# Example usage:
# Assuming taxident.all is your data.table and ranks is the vector of columns to concatenate
# paste_columns_dt(taxident.all, ranks, '::', 'lineage', 1)

# Example usage:
# paste_columns(dt, c("col1", "col2", "col3"), "-", "new_pasted_col", 2)
