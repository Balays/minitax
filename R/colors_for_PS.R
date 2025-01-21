

## get colors for given ps object

colors_for_PS <- function(ps, rank, pal.man, unassigned='missing', unknowns='unknown_') {

  if(rank != 'taxa_names') {
    colors   <- setdiff(as.character(ps@tax_table[, rank]), unassigned)
  } else {
    colors   <- setdiff(as.character(taxa_names(ps)), unassigned)
  }

  unknown_cols <- grep(unknowns, colors, value = T)
  unknown_cols <- as.character(grey.colors(length(unknown_cols)))
  names(unknown_cols) <- grep(unknowns, colors, value = T)

  unknown_cols <- c(unknown_cols, unassigned = 'black')

  colors <- grep(unassigned, colors, value = T, invert = T)

  colors <- grep(unknowns, colors, value = T, invert = T)

  my.pal <- pal.man[1:length(colors)]

  names(my.pal) <- colors

  my.pal <- c(my.pal, unknown_cols)

}


#my.pal <- cols_for_ps(ps=ps.glom.norm.top.merged, rank=t, pal.man)
