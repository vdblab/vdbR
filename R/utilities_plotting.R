# Currently this is just hard coded, ideally would read in from postgres:
base_colors = c(other="#BBBBBB",
                p__Actinobacteria = "#A77097", p__Bacteroidetes = "#51AB9B", p__Proteobacteria = "#AB3535", # "#FF0000". Ying had used full red spectrum. I used DD0000 to increase spectrum.
                o__Clostridiales = "#9C854E",
                o__Erysipelotrichales = "#7A6920", #o__Erysipelotrichales = "#ACA54D", #Erysipelotrichales order has been added mainly due to high frequency of genus: Erysipelatoclostridium.
                f__Lachnospiraceae = "#EC9B96", f__Ruminococcaceae = "#FF7106" ,#"#9AAE73",
                g__Enterococcus = "#129246", g__Lactobacillus = "#3B51A3", g__Staphylococcus = "#F1EB25", g__Streptococcus = "#9FB846",
                g__Escherichia = "#800026", g__Klebsiella = "#FF2211",
                g__Akkermansia = "#377EB8");

taxa_order = c(kingdom = "k", phylum="p", class="c", order="o", family="f", genus="g")

#' This function will generate a color close to the original color comp provided
#' 
#' @param rgb_component the component to randomly adjust
#' @param magnitude the maximum amount to adjust rgb component by (default 4)
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



#' This function will generate colors for the taxa_levels provided based on the base_colors and the taxa levels.
#' New colors will be generated in a spread on the appropriate rank. 
#' 
#' @param palette the palette to overwrite colors for.
#' @param full_taxonomy full taxonomy rows for the species to generate colors for.
#' @param rank the rank of the palette
#' @export
#' @name rename_taxa_colors
rename_taxa_colors <- function(palette, full_taxonomy, rank){
  color_group_tax_info = data.frame(group = names(base_colors)) %>%
    dplyr::mutate(tax_level = gsub("(.*)__.*", "\\1", group))
  inds_to_adjust = rep(FALSE, length(palette))
  new_palette = unname(palette)
  shuffled_names = rownames(full_taxonomy)
  
  # For all the taxa names, update the color and "blur" the colors from higher clades.
  for(tax_level in colnames(full_taxonomy) ){
    new_palette[inds_to_adjust] <- lapply(new_palette[inds_to_adjust], blur_color)
    for(clade in unique(full_taxonomy[tax_level][[1]])){
      if (clade %in% names(base_colors)){
        new_palette[which(full_taxonomy[tax_level] == clade)] <- base_colors[clade]
        inds_to_adjust = inds_to_adjust | full_taxonomy[tax_level] == clade
        s_index <- grepl(clade, shuffled_names)
        shuffled_names <- c(shuffled_names[s_index], shuffled_names[!s_index])
      }
    }
  }
  
  #This will swap around the order of the names so that groups of related genera are together in the palette:
  ordered_palette <- data.frame(
    color = list_c(new_palette),
    name = full_taxonomy[rank][[1]],
    full_tax = rownames(full_taxonomy)
  ) %>% left_join(data.frame(
    full_tax = shuffled_names,
    order_index = seq(1, length(palette))
  )) %>%
    arrange(order_index)
  
  # Make the new palette and add in the "other" color:
  new_palette <- ordered_palette$color
  new_palette <- c(new_palette, "lightgrey")
  new_palette <- as.character(new_palette)
  
  names(new_palette) <- c(ordered_palette$name, "Other")
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
#' @export
#' @name make_microviz_palette
make_microviz_palette <- function(phy_seq_obj, n, rank){
  
  # for our list of required packages run through and insure they are installed:
  req_packages = c("microViz")
  for (package in req_packages){
    if (!requireNamespace(package, quietly=TRUE)){
      warning(paste("The", package, "package must be installed in order to make use of the make_microviz_palette function."))
      return(NULL)
    }
  }
  
  #make the default palette based on phyloseq object and params:
  my_palette <- microViz::tax_palette(phy_seq_obj, n=n, rank=rank, add = NA) #get a palette with n/rank
  index_rank <- grep(rank, colnames(phyloseq::tax_table(phy_seq_obj)))
  
  taxonomy <- phy_seq_obj %>%
    phyloseq::tax_table() %>%
    as.data.frame() %>%
    dplyr::select(1:all_of(index_rank)) %>%
    dplyr::filter(.[[index_rank]] %in% names(my_palette)) %>%
    unique()
  
  #rename colors using basenames:
  return(rename_taxa_colors(my_palette, taxonomy, rank))
  
}
