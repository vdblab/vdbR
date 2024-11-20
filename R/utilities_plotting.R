#' This function will create a microviz palette for phy_seq_obj and use db of colors where possible. 
#' Note, this function will use Genus names as the keys for the color table. If working with SGB names
#' it is recommended to first clean the tax table of the phyloseq object with the clean_SGB_genus 
#' function provided in this package.
#'
#' @param phy_seq_obj the phyloseq object for which the palette is being made. 
#' @param n the number of unique taxa to use in the palette (should match what is being used in the plot)
#' @param rank the taxonomic rank being used in the plot and the palette 
#' @param color_db a table to use for overwriting color values. requires genus, color and color_label_group as columns. (eg. asv_annotation_blast_color_ag)
#' @export
#' @name make_microviz_palette
#' 
make_microviz_palette <- function(phy_seq_obj, n, rank, color_db){
  
  # for our list of required packages run through and insure they are installed:
  req_packages = c("microViz")
  for (package in req_packages){
    if (!requireNamespace(package, quietly=TRUE)){
      warning(paste("The", package, "package must be installed in order to make use of the make_microviz_palette function."))
      return(NULL)
    }
  }
  
  # test that the genus column is formatted consistently for color_db
  test_df <- phyloseq::tax_table(phy_seq_obj) %>%
    as.data.frame() %>%
    dplyr::mutate(formatted = grepl("^[gfocp]__.+", genus) | grepl("other", genus))
  
  if(!all(test_df$formatted)){
    warning("Not all entries in provided genus column of phy_seq_obj seem to be \
            Properly formatted,  Consider cleaning using clean_SGB_genus. \
            Failure to do so means colors may not be successfully overwritten. ")
  }
  
  # process color db, making it unique with respect to genus and selecting required columns:
  color_table <- color_db %>%
    as.data.frame() %>%
    dplyr::group_by(genus) %>%
    dplyr::mutate(first_seen = dplyr::row_number() == 1) %>%
    dplyr::ungroup() 
  color_table <- color_table[which(color_table$first_seen),] %>% #drop all non-unique genus entries. 
    dplyr::select(genus, color_label_group, color)
  
  #make the default palette based on phyloseq object and params:
  my_palette <- microViz::tax_palette(phy_seq_obj, n=n, rank=rank) #get a palette with n/rank
  
  #convert palette names into color labels.  Note we strip any numbers here, and remove all unclassified tags.
  color_labels <- names(my_palette) %>% 
    as.data.frame()
  names(color_labels)[1] <- "original_name"
  color_labels <- color_labels %>%
    dplyr::mutate(color_label_group = lapply(strsplit(original_name, ' '), '[', 1)) %>%
    dplyr::mutate(color_label_group = stringr::str_remove(color_label_group, '_unclassified'))
  
  #for values in the palette which have another color in the db, overwrite:
  for (i in 1:nrow(color_labels)){
    genus = ''
    name = color_labels[i,]$color_label_group
    if (startsWith(name, 'g__')){ 
      genus <- substr(name, 4, nchar(name))
    }# prefer overwriting color by the genus name:
    if (genus != '' & genus %in% color_table$genus) {
      new_color = color_table[which(color_table$genus == genus),]$color
      my_palette[color_labels[i,]$original_name] = new_color
    } # when that is not possible, overwrite by first color_label_group:
    else if (color_labels[i,]$color_label_group %in% color_table$color_label_group){
      new_color = color_table[which(color_table$color_label_group == color_labels[i,]$color_label_group),][1,]$color
      my_palette[color_labels[i,]$original_name] = new_color
    }
  }
  return(my_palette)
}
