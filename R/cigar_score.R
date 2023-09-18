require(purrr)


cigar_score <- function(cigar, match_score    =  1, 
                        mismatch_score        = -3,
                        insertion_score       = -2, 
                        deletion_score        = -2, 
                        gap_opening_penalty   = -2, 
                        gap_extension_penalty = -1) {
  
  # Helper function to calculate score for each operation
  calc_score <- function(op) {
    num_bases <- as.integer(gsub("[MIDS]", "", op))
    op_type <- substr(op, nchar(op), nchar(op))
    
    switch(op_type,
           M = match_score * num_bases,
           I = insertion_score + gap_opening_penalty + gap_extension_penalty * (num_bases - 1),
           D = deletion_score + gap_opening_penalty + gap_extension_penalty * (num_bases - 1),
           S = 0,
           mismatch_score * num_bases
    )
  }
  
  # Helper function to calculate alignment length for each operation
  calc_length <- function(op) {
    as.integer(gsub("[MIDS]", "", op))
  }
  
  # Split the CIGAR string into separate operations
  operations <- strsplit(gsub("([MIDS])", "\\1 ", cigar), " ")[[1]]
  operations <- operations[operations != ""]
  
  # Calculate the score and alignment length using map_dbl
  score <- sum(map_dbl(operations, calc_score))
  alignment_length <- sum(map_dbl(operations, calc_length))
  
  # Normalize by alignment length
  normalized_score <- score / alignment_length
  
  return(normalized_score)
}
