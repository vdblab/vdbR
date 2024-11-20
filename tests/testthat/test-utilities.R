warning_about_singlecol_metadata <-"Phyloseq object must have at least one metadata column in addition to the sampleids; creating a new columns from the sampleids called x"
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

test_that("making a phyloseq object from some metadata works", {
  ps <- vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
  expect_equal(phyloseq::rank_names(ps), c("kingdom", "phylum", "class", "ordr", "family", "genus", "species"))
})
test_that("making a phyloseq object without a metadata dataframe fails", {
  # this is a vector of sampleids
  expect_error(
    vdb_make_phylo(test_metadata[, c("sampleid")], sampleid_col = "sampleid")
  )
})
test_that("making a phyloseq object from a single metadata column dataframe by adding dummy column", {
  ps <- suppressWarnings(vdb_make_phylo(tibble::as_tibble(test_metadata)[, c("sampleid")], sampleid_col = "sampleid"))
  expect_equal(phyloseq::rank_names(ps), c("kingdom", "phylum", "class", "ordr", "family", "genus", "species"))
  expect_equal(phyloseq::sample_data(ps)$x, phyloseq::sample_names(ps))
})


test_that("check handling metaphlan sample names single sequencing run", {
  data = data.frame(
    sample_id = c("PJlibI_MR8R", "PJlibI_MR7R")
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=92), warning_about_singlecol_metadata)
  expect_equal(sort(phyloseq::sample_names(ps)), sort(c("15423_PJlibI_MR8R_H7NNYDSX7", "15423_PJlibI_MR7R_H7NNYDSX7")))
})

test_that("check handling metaphlan sample names mixed sequencing run", {
  data = data.frame(
    sample_id = c("3384K", "3384A")
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=66), warning_about_singlecol_metadata)
  expect_equal(sort(phyloseq::sample_names(ps)), sort(data$sample_id))
})

test_that("check handling metaphlan sample names just multiple sequencing run", {
  data = data.frame(
    sample_id = c("3384K", "3384J")
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=66), warning_about_singlecol_metadata)
  expect_equal(sort(phyloseq::sample_names(ps)), sort(data$sample_id))
})

test_that("process_metadata doesn't allow NAs", {
  data = data.frame(
    sample_id = c("3384K", "3384J", NA)
  )
  
  expect_error(process_metadata(data, sampleid_col = "sample_id"))
})
test_that("process_metadata doesn't allow duplicate sample IDs", {
  data = data.frame(
    sample_id = c("3384K", "3384J", "3384J")
  )
  expect_error(process_metadata(data, sampleid_col = "sample_id"))
})
test_that("The microviz palette making function correctly overwrites color.", {
  if (requireNamespace('microViz', quietly=TRUE)){ #only run test if microViz is installed.
    ps <- vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
    phyloseq::tax_table(ps) <- microViz::tax_fix(ps) 
    #rename genus to match the keys in the color db
    phyloseq::tax_table(ps) <- phyloseq::tax_table(ps) %>%
      as.data.frame() %>%
      dplyr::mutate(genus = ifelse(grepl("^[gfocp]__.+", genus), genus, paste0("g__", genus))) %>%
      as.matrix() %>%
      phyloseq::tax_table()
    default_palette = microViz::tax_palette(ps, n=15, rank='genus')
    connect_database()
    get_table_from_database('asv_annotation_blast_color_ag')
    adj_palette = make_microviz_palette(ps, n=15, rank='genus', color_db=asv_annotation_blast_color_ag)
    expect(as.character(adj_palette["g__Enterococcus"]) != as.character(default_palette['g__Enterococcus']),
           failure_message="Enterococcus record was not correctly overwritten.")
  }
})
test_that("test that clean_SGB_names successfully overwrites SGB names at various levels", {
  taxmat = matrix(sample(letters, 35, replace = TRUE), nrow = 5, ncol = 7)
  taxmat[,2] = sample(letters, 5, replace = FALSE) #make sure we are unique at Phylum
  taxmat[2:5,6] = 'GGB'
  taxmat[3:5,5] = 'FGB'
  taxmat[4:5,4] = 'OFG'
  taxmat[5,3] = 'CFG'
  rownames(taxmat) <- paste0("OTU", 1:nrow(taxmat))
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  clean_mat <- clean_SGB_genus(phyloseq::tax_table(taxmat))
  expect(startsWith(clean_mat[1,6], 'g__'), failure_message = "defaulting to genus not working")
  expect(startsWith(clean_mat[2,6], 'f__'), failure_message = "family overwrite not working")
  expect(startsWith(clean_mat[3,6], 'o__'), failure_message = "order overwrite not working")
  expect(startsWith(clean_mat[4,6], 'c__'), failure_message = "class overwrite not working")
  expect(startsWith(clean_mat[5,6], 'p__'), failure_message = "phylum overwrite not working")
})
test_that("test clean_SGB_genus name reformatting", {
  taxmat = matrix(sample(letters, 35, replace = TRUE), nrow = 5, ncol = 7)
  rownames(taxmat) <- paste0("OTU", 1:nrow(taxmat))
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  clean_mat <- clean_SGB_genus(phyloseq::tax_table(taxmat))
  lowercase_names <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  cleaned_names = colnames(phyloseq::tax_table(clean_mat))
  expect(setequal(intersect(lowercase_names, cleaned_names), cleaned_names),
         failure_message = "The expected set of cleaned up names is not a subset of columns returned. ")
})


