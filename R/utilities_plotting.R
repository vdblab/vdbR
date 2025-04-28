taxa_order = c(kingdom = "k", phylum="p", class="c", order="o", family="f", genus="g")

#' This function will generate a color close to the original color comp provided
#' 
#' @param rgb_component the component to randomly adjust
#' @param magnitude the maximum amount to adjust rgb component by (default 15)
#' @export
#' @name norm_generate_component
norm_generate_component <- function(rgb_component, magnitude = 15){
  rgb_component <- rgb_component + sample(-magnitude:magnitude, 1)
  rgb_component <- min(255, max(0, rgb_component))
}

#' This function will generate a color close to the original color provided
#' @importFrom grDevices col2rgb rgb
#' @param color color to randomly adjust (as hex color)
#' @export
#' @name blur_color
blur_color <- function(color){
  rgb_comps <- col2rgb(color)
  return(rgb(
    norm_generate_component(rgb_comps[1]), 
    norm_generate_component(rgb_comps[2]), 
    norm_generate_component(rgb_comps[3]), 
    maxColorValue = 255, 
    names = NULL))
  
}



#' This function will generate colors for the taxa_levels provided based on the base_palette provided and the taxa levels.
#' New colors will be generated in a spread on the appropriate rank. 
#' 
#' @param palette the palette to overwrite colors for.
#' @param full_taxonomy full taxonomy rows for the species to generate colors for.
#' @param rank the rank of the palette
#' @param base_palette the base palette to use to generate taxonomy specific colors (named)
#' @export
#' @name rename_taxa_colors
rename_taxa_colors <- function(palette, full_taxonomy, rank, base_palette, shuf_genus=T){
  color_group_tax_info = data.frame(group = names(base_palette)) %>%
    dplyr::mutate(tax_level = gsub("(.*)__.*", "\\1", group))
  inds_to_adjust = rep(FALSE, length(palette))
  new_palette = unname(palette)
  shuffled_names <- full_taxonomy %>%
    dplyr::mutate(joined = paste(phylum, class, order, family, genus)) %>%
    dplyr::pull(joined)
  
  # For all the taxa names, update the color and "blur" the colors from higher clades.
  for(tax_level in colnames(full_taxonomy) ){
    new_palette[inds_to_adjust] <- lapply(new_palette[inds_to_adjust], blur_color)
    for(clade in unique(full_taxonomy[tax_level][[1]])){
      if (clade %in% names(base_palette)){
        new_palette[which(full_taxonomy[tax_level] == clade)] <- base_palette[clade]
        inds_to_adjust = inds_to_adjust | full_taxonomy[tax_level] == clade
        s_index <- grepl(clade, shuffled_names)
        if (shuf_genus | !tax_level == "genus"){ # don't shuffle genus level is the flag is true to not mess up order
          shuffled_names <- c(shuffled_names[s_index], shuffled_names[!s_index])
        }
      }
    }
  }
  
  
  #This will swap around the order of the names so that groups of related genera are together in the palette:
  ordered_palette <- data.frame(
    color = unlist(new_palette),
    name = full_taxonomy[rank][[1]],
    full_tax = rownames(full_taxonomy) 
  ) %>% 
    dplyr::mutate(joined = paste(full_taxonomy$phylum, full_taxonomy$class, full_taxonomy$order, full_taxonomy$family, full_taxonomy$genus))  %>% 
    dplyr::left_join(data.frame(
    joined = shuffled_names,
    order_index = seq(1, length(palette))
  )) %>% 
    dplyr::arrange(order_index)
  
  # Make the new palette:
  new_palette <- ordered_palette$color
  new_palette <- as.character(new_palette)
  names(new_palette) <- ordered_palette$name
  
  gray_locs <- grepl("gray", new_palette)
  new_palette <- c(new_palette[!gray_locs], new_palette[gray_locs])
  new_palette <- c(new_palette, c(Other = "gray"))
  
  return(new_palette)
}

#' This function will create a microviz palette for phy_seq_obj and use db of colors where possible. 
#' Note, this function will use Genus names as the keys for the color table. If working with SGB names
#' it is recommended to first clean the tax table of the phyloseq object with the clean_SGB_genus 
#' function provided in this package.
#'
#' @param phy_seq_obj the phyloseq object for which the palette is being made. 
#' @param n the number of unique taxa to use in the palette (should match what is being used in the plot)
#' @param rank the taxonomic rank being used in the plot and the palette 
#' @param taxo_palette (optional) the base palette to use to generate taxonomy specific colors (named). Defaults to Ying Taur palette
#' @export
#' @name make_microviz_palette
make_microviz_palette <- function(phy_seq_obj, n, rank, taxo_palette=NA, shuf_genus=T){
  
  # for our list of required packages run through and insure they are installed:
  req_packages = c("microViz")
  for (package in req_packages){
    if (!requireNamespace(package, quietly=TRUE)){
      warning(paste("The", package, "package must be installed in order to make use of the make_microviz_palette function."))
      return(NULL)
    }
  }
  
  #make a default palette of the appropriate size with default color gray:
  top_taxa <- microViz::tax_top(phy_seq_obj, n=n, rank=rank)
  my_palette <- rep('gray', n)
  names(my_palette) <- c(top_taxa)
  index_rank <- grep(rank, colnames(phyloseq::tax_table(phy_seq_obj)))
  
  taxonomy <- phy_seq_obj %>%
    phyloseq::tax_table() %>%
    as.data.frame() %>%
    dplyr::select(1:dplyr::all_of(index_rank)) %>%
    dplyr::filter(.[[index_rank]] %in% names(my_palette)) %>%
    unique()
  
  if (!is.character(taxo_palette)){
    #if no palette is provided will default to the Ying Taur palette.
    taxo_palette = vdbR::full_palette
  }
  
  #rename colors using basenames:
  return(rename_taxa_colors(my_palette, taxonomy, rank, taxo_palette, shuf_genus))
  
}
