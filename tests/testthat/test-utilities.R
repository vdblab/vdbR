warning_about_singlecol_metadata <-"Phyloseq object must have at least one metadata column in addition to the sampleids; creating a new columns from the sampleids called x"
error_about_only_supporting_biobakery <- "Currently only Biobakery/Metaphlan based make_phylo_mgx is supported. Please open a ticket on Github if you need extended functionality! "
error_about_sampleid_col_not_in_metadata <- "sampleid_col not in colnames of metadata supplied."
error_about_mixed_sample_and_experiment_level_analyses <- paste0("The same experiments are in multiple analyses. ", 
        "This can happen when someone runs an app both per-sample and per-experiment. ",
        "Resolve this in Isabl to avoid confusion, or rerun this code with ",
        "choose_max_experiment = True")
test_metadata <- structure(
  list(
    sampleid=c("1143N", "1252XX", "1303G", "1379N", "142E", "1528W", "1562G", 
               "1567C", "224C", "430N", "725E", "729I", "933A", "966E", "FMT.0011C", 
               "FMT.0012P", "FMT.0073D", "FMT.0144C", "FMT.0212D", "FMT.0241C"
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
  skip_if(Sys.getenv("GITHUB_ACTIONS") != "")
  connect_database(bundled = FALSE)
  ps <- vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
  expect_equal(phyloseq::rank_names(ps), c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
})

test_that("making a phyloseq object from some metadata works from testdb", {
  connect_database(bundled = TRUE)
  ps <- vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
  expect_equal(phyloseq::rank_names(ps), c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
})



test_that("making a phyloseq object without a metadata dataframe fails", {
  # this is a vector of sampleids
  expect_error(
    vdb_make_phylo(test_metadata[, c("sampleid")], sampleid_col = "sampleid")
  )
})

test_that("making a phyloseq object from a single metadata column dataframe by adding dummy column", {
  connect_database(bundled = TRUE)
  ps <- suppressWarnings(vdb_make_phylo(tibble::as_tibble(test_metadata)[, c("sampleid")], sampleid_col = "sampleid"))
  #expect_equal(phyloseq::rank_names(ps), c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
  expect_equal(phyloseq::sample_data(ps)$x, phyloseq::sample_names(ps))
})

# app_id = 1 = Fake Biobakery application in test_tables:
test_that("running vdb_make_phylo_mgx fails if run on non-biobakery app", {
  data = data.frame(
    sample_id = c("PJlibI_MR8R", "PJlibI_MR7R")
  )
  expect_error(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=2, testing = TRUE), error_about_only_supporting_biobakery)
})

test_that("vdb_make_phylo_mgx fails sampleid_col passed is not in the metadata", {
  data = data.frame(
    sample_id = c("PJlibI_MR8R", "PJlibI_MR7R")
  )
  expect_error(ps <- vdb_make_phylo_mgx(data, sampleid_col = "something_wrong_here_not_valid_name", app_id=1, testing = TRUE),
    error_about_sampleid_col_not_in_metadata)
})

test_that("check handling metaphlan sample names single sequencing run", {
  data = data.frame(
    sample_id = c("simple_sample_1_exp")
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=1, testing = TRUE), warning_about_singlecol_metadata)
  expect_equal(sort(phyloseq::sample_names(ps)), sort(data$sample_id))
})

test_that("check handling metaphlan sample names mixed sequencing run (many experiments run sample level)", {
  data = data.frame(
    sample_id = c("sample_many_exps_run_correctly")
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=1, testing = TRUE), warning_about_singlecol_metadata)
  expect_equal(sort(phyloseq::sample_names(ps)), sort(data$sample_id))
})

test_that("check handling metaphlan both single experiment and multiple experiment run correctly", {
  data = data.frame(
    sample_id = c("sample_many_exps_run_correctly", "simple_sample_1_exp")
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=1, testing = TRUE), warning_about_singlecol_metadata)
  expect_equal(sort(phyloseq::sample_names(ps)), sort(data$sample_id))
})

test_that("error if run on sample with mixed experiment level and sample level results", {
  data = data.frame(
    sample_id = c("sample_many_exps_mixed_analyses")
  )
  expect_error(expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=1, testing = TRUE), warning_about_singlecol_metadata), error_about_mixed_sample_and_experiment_level_analyses)
})

test_that("mixed sample_level and experiment level analyses will work if choose_max_experiment = TRUE", {
  data = data.frame(
    sample_id = c("sample_many_exps_mixed_analyses")
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=1, choose_max_experiment = TRUE, testing = TRUE), warning_about_singlecol_metadata)
  expect_equal(sort(phyloseq::sample_names(ps)), sort(data$sample_id))
})


test_that("real data mixed sample_level and experiment level analyses will work if choose_max_experiment = TRUE", {
  data = data.frame(
    sample_id = c("PJlibH_GutZymo", "GutZymoDiet", "GUTZYMO.1177", "MMF_hc_E_GutZymo", 
                  "QM154_GutZymo", "ZB132_GutZymo", "PF53_GutZymo", "PC26_GutZymo", 
                  "GutZymo", "GutZymo_28", "Zymo", "zymostandard2", "zymostandard1"
    )
  )
  expect_warning(ps <- vdb_make_phylo_mgx(data, sampleid_col = "sample_id", app_id=66, choose_max_experiment = TRUE, testing = FALSE), warning_about_singlecol_metadata)
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


