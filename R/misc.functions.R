
dup <- function (x) x[duplicated(x)]

luniq <- function (x) length(unique(x))

is.0  <- function(x) x == 0

## add prime 3s and prime 5s based on start, end and strand
add.primes <- function(df) {
    
  
  df$prime5[df$strand == '+' & !is.na(df$strand)] <- df$start[df$strand == '+' & !is.na(df$strand) ]
  df$prime3[df$strand == '+' & !is.na(df$strand)] <- df$end[df$strand   == '+' & !is.na(df$strand) ]
  ##
  df$prime5[df$strand == '-' & !is.na(df$strand)] <- df$end[df$strand   == '-' & !is.na(df$strand) ]
  df$prime3[df$strand == '-' & !is.na(df$strand)] <- df$start[df$strand == '-' & !is.na(df$strand) ]
  
  return(df)
  
}


write_tsv <- function(x, filename, ...) {
  
  write.table(x, filename, row.names = F, quote = F, sep='\t', ...)
  
}



Mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

