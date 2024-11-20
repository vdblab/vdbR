warning_about_singlecol_metadata <-"Phyloseq object must have at least one metadata column in addition to the sampleids; creating a new columns from the sampleids called x"
error_about_only_supporting_biobakery <- "Currently only Biobakery/Metaphlan based make_phylo_mgx is supported. Please open a ticket on Github if you need extended functionality! "
error_about_sampleid_col_not_in_metadata <- "sampleid_col not in colnames of metadata supplied."
error_about_mixed_sample_and_experiment_level_analyses <- paste0("The same experiments are in multiple analyses. ", 
        "This can happen when someone runs an app both per-sample and per-experiment. ",
        "Resolve this in Isabl to avoid confusion, or rerun this code with ",
        "choose_max_experiment = True")
test_metadata <- structure(
  list(
    sampleid = c(
      "2133C", "1773A", "1773I", "90.tp.51",
      "223D", "2065C", "1773O", "1773M", "2133D", "2065D", "90.tp.97",
      "1773E", "2065B", "2065E", "90.tp.57", "2133G", "2065F", "1773G",
      "FMT.0113I", "2065G"
    ),
    Conditioning = c(
      "Flu/BU4", "Flu/BU4",
      "Flu/BU4", "Flu/BU3", "Flu/BU4", "Flu/BU4", "Flu/BU4", "Flu/BU4",
      "Flu/BU4", "Flu/BU4", "Flu/BU4", "Flu/BU4", "Flu/BU4", "Flu/BU4",
      "Flu/BU4", "Flu/BU4", "Flu/BU4", "Flu/BU4", "Flu/BU4", "Flu/BU4"
    ),
    source = c(
      "PBSC", "PBSC", "PBSC", "PBSC", "PBSC", "PBSC",
      "PBSC", "PBSC", "PBSC", "PBSC", "PBSC", "PBSC", "PBSC", "PBSC",
      "PBSC", "PBSC", "PBSC", "PBSC", "PBSC", "PBSC"
    ),
    age = c(
      "60 and over ",
      "60 and over ", "60 and over ", "60 and over ", "60 and over ",
      "40-59", "60 and over ", "60 and over ", "60 and over ", "40-59",
      "40-59", "60 and over ", "40-59", "40-59", "40-59", "60 and over ",
      "40-59", "60 and over ", "40-59", "40-59"
    ),
    fakemrn = c(
      "00000001",
      "00000002", "00000003", "00000004", "00000005", "00000001", "00000002",
      "00000003", "00000004", "00000005", "00000001", "00000002", "00000003",
      "00000004", "00000005", "00000001", "00000002", "00000003", "00000004",
      "00000005"
    )
  ),
  row.names = c(NA, -20L), class = "data.frame"
) %>% 
  dplyr::slice(1:5)


test_that("The microviz palette making function correctly overwrites color.", {
  #only run test if microViz is installed. and on-prep
  skip_if(Sys.getenv("GITHUB_ACTION") != "")
  if (requireNamespace('microViz', quietly=TRUE) ){ 
    connect_database()
    ps <- vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
    phyloseq::tax_table(ps) <- microViz::tax_fix(ps) 
    #rename genus to match the keys in the color db
    phyloseq::tax_table(ps) <- phyloseq::tax_table(ps) %>%
      as.data.frame() %>%
      dplyr::mutate(genus = ifelse(grepl("^[gfocp]__.+", genus), genus, paste0("g__", genus))) %>%
      as.matrix() %>%
      phyloseq::tax_table()
    default_palette = microViz::tax_palette(ps, n=15, rank='genus')
    adj_palette = make_microviz_palette(ps, n=15, rank='genus')
    expect(as.character(adj_palette["g__Enterococcus"]) != as.character(default_palette['g__Enterococcus']),
           failure_message="Enterococcus record was not correctly overwritten.")
  }
})


