cigar_to_length <- function(cigar) {
  operations <- unlist(str_extract_all(cigar, "[MIDNSHP=X]"))
  lengths    <- as.numeric(unlist(str_extract_all(cigar, "\\d+")))
  
  consuming_ops <- c("M", "I", "S", "=", "X")
  sum(lengths[operations %in% consuming_ops])
}