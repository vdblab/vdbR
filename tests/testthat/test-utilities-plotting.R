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


test_that("The microviz palette making function correctly changes genera that should be changed.", {
  #only run test if microViz is installed. and on-prem
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


test_that("The microviz palette looks similar to previous palette.", {
  #only run test if on-prem
  skip_if(Sys.getenv("GITHUB_ACTION") != "")

  if (requireNamespace('microViz', quietly=TRUE) ){ 
    connect_database()
    data = data.frame(
      sample_id = c("sample_many_exps_mixed_analyses")
    )
    expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=1, choose_max_experiment = TRUE, testing = TRUE) %>%
      microViz::tax_fix() %>%
      phyloseq::transform_sample_counts(function(x) round(x*1000000, 0) ))
    
    phyloseq::tax_table(ps) <- clean_SGB_genus(ps)
    n_taxa = 15
    rank = "genus"
    
    set.seed(123)
    test_palette <- make_microviz_palette(ps, n_taxa, rank, shuf_genus = F)
    expected_palette <- readRDS("palette_for_testthat.rds")
    expect_equal(names(test_palette), names(expected_palette))
    genera = names(test_palette)
    expect_equal(test_palette[genera], expected_palette[genera])
  }
})

test_that("The palette is overwritten when custom palette is passed.", {
  #only run test if on-prem
  skip_if(Sys.getenv("GITHUB_ACTION") != "")
  if (requireNamespace('microViz', quietly=TRUE) ){ 
    connect_database()
    data = data.frame(
      sample_id = c("sample_many_exps_mixed_analyses")
    )
    expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=1, choose_max_experiment = TRUE, testing = TRUE) %>%
                     microViz::tax_fix() %>%
                     phyloseq::transform_sample_counts(function(x) round(x*1000000, 0) ))
    
    phyloseq::tax_table(ps) <- clean_SGB_genus(ps)
    n_taxa = 15
    rank = "genus"
    
    set.seed(123)
    test_palette <- make_microviz_palette(ps, n_taxa, rank, c(g__Streptococcus = "blue"))
    expect_equal(test_palette["g__Streptococcus"], c(g__Streptococcus ="blue"))
    expect_equal(test_palette["g__Bacteroides"], c(g__Bacteroides ="gray"))
  }
})


test_that("sample legacy composition plots looks similiar", {
  #only run test if on-prem
  skip_if(Sys.getenv("GITHUB_ACTION") != "")
  connect_database()
    get_table_from_database("frozen_set_ag")
  get_table_from_database("asv_annotation_blast_color_ag")
  set.seed(123)
  phy_amplicon <- vdb_make_phylo(frozen_set_ag %>% dplyr::sample_n(5)) %>% 
    phyloseq::tax_glom("genus")  %>% 
    microViz::tax_fix()
  
  phyloseq::tax_table(phy_amplicon) <- clean_refseq_genus(phy_amplicon)
  # get sample order  based on the higherarchical clustering of the bray curtis distances
  sample_order <- labels(as.dendrogram(hclust(phyloseq::distance(phy_amplicon, method = "bray"))))
  fullpal <- make_microviz_palette(phy_amplicon, n=phyloseq::ntaxa(phy_amplicon), rank="genus", shuf_genus = FALSE) 
  
  p_comp <- phy_amplicon %>%
    microViz::comp_barplot(tax_level = "genus", palette = fullpal, n_taxa=phyloseq::ntaxa(phy_amplicon), bar_outline_colour = NA, sample_order = sample_order, tax_order = names(fullpal) ) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) + 
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal") + 
    ggplot2::guides(fill="none") + 
    ggplot2::scale_y_continuous(expand=c(0, 0)) 
  
  # classic plot
  # merge with color palete
  bar_data <- phy_amplicon %>% phyloseq::psmelt() %>% 
    dplyr::left_join(asv_annotation_blast_color_ag %>% dplyr::select(OTU=asv_key, color)) 
  # plot after relative abundance conversion and genus-level aggregation
  p_old_comp <- ggplot2::ggplot(
    bar_data %>% 
      dplyr::group_by(Sample) %>%
      dplyr::mutate(Abundance=Abundance/sum(Abundance)) %>% 
      dplyr::group_by(Sample, color, family,genus) %>%  
      dplyr::summarize(Abundance=sum(Abundance)), ggplot2::aes(x=Sample, y=Abundance, fill=color)) +
    ggplot2::geom_bar(stat="identity", color=NA, width = 1) + 
    ggplot2::scale_fill_identity(guide="none")  +
    ggplot2::scale_x_discrete(limits=sample_order) + 
    ggplot2::scale_y_continuous(expand=c(0,0))  + 
    ggplot2::labs(title="old") 
  pdf(file.path( "old_vs_new_comp_barplots.pdf"), width = 8, height = 10)
  print(p_comp  + ggplot2::labs(title="new") )
  print(p_old_comp)
  dev.off()
testthat::expect_snapshot_file("old_vs_new_comp_barplots.pdf")
testthat::expect_message(message("please check the temporary plots for equivalency in tests/testthat/utilities-plotting/old_vs_new_comp_barplots.pdf "))
})
