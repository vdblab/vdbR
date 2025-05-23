% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities_plotting.R
\name{rename_taxa_colors}
\alias{rename_taxa_colors}
\title{This function will generate colors for the taxa_levels provided based on the base_palette provided and the taxa levels.
New colors will be generated in a spread on the appropriate rank.}
\usage{
rename_taxa_colors(palette, full_taxonomy, rank, base_palette, shuf_genus = T)
}
\arguments{
\item{palette}{the palette to overwrite colors for.}

\item{full_taxonomy}{full taxonomy rows for the species to generate colors for.}

\item{rank}{the rank of the palette}

\item{base_palette}{the base palette to use to generate taxonomy specific colors (named)}

\item{shuf_genus}{(Boolean, optional) if false will not shuffle order of colors at genus level (useful when using a hierarchical palette with many colors redefined at the genus level). Default T.}
}
\description{
This function will generate colors for the taxa_levels provided based on the base_palette provided and the taxa levels.
New colors will be generated in a spread on the appropriate rank.
}
