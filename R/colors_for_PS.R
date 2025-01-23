colors_for_PS <- function(ps=ps.glom.norm.top.merged, rank='genus', pal.man, unassigned='missing', unknowns='unknown_') {

  if(rank != 'taxa_names') {
    colors   <- setdiff(as.character(ps@tax_table[, rank]), unassigned)
  } else {
    colors   <- setdiff(as.character(taxa_names(ps)), unassigned)
  }

  unknown_cols <- grep(unknowns, colors, value = T)
  unknown_cols <- as.character(grey.colors(length(unknown_cols)))
  names(unknown_cols) <- grep(unknowns, colors, value = T)

  unknown_cols <- c(unassigned = 'black', unknown_cols)
  names(unknown_cols)[1] <- unassigned

  colors <- grep(unassigned, colors, value = T, invert = T)

  colors <- grep(unknowns, colors, value = T, invert = T)

  my.pal <- pal.man[1:length(colors)]

  names(my.pal) <- colors

  my.pal <- c(my.pal, unknown_cols)

}
