#' Ensure that metadata is a dataframe without NA sampleid or duplicated sampleids.
#' cast metadata to dataframe so we can do convenient indexing 
#'
#' @param metadata metadata dataframe
#' @param sampleid_col name of column in metadata containing sample_ids 
#' @name process_metadata
#'
process_metadata <- function(metadata, sampleid_col){
  metadata <- as.data.frame(metadata)
  if(length(unique(metadata[, sampleid_col])) != nrow(metadata)){
    stop(paste0("Metadata's sample id column ", sampleid_col, " must have unique entries"))
  } else if(any(is.na(metadata[, sampleid_col])) ){
    stop(paste0("Metadata's sample id column ", sampleid_col, " must not have NAs"))
  }
  if (ncol(metadata) == 1){
    warning("Phyloseq object must have at least one metadata column in addition to the sampleids; creating a new columns from the sampleids called x")
    metadata$x <- metadata[, sampleid_col]
  }
  metadata
}


#A little helper function which just renames taxonomic names to be consistent:
convert_taxo_names_to_lowercase <- function(x){
  lowercase_names <- c('class', 'domain', 'family', 'genus', 'order', 'phylum', 'species')
  for (col in colnames(x)){
    if (tolower(col) %in% lowercase_names & !(col %in% lowercase_names)){
      names(x)[names(x) == col] <- tolower(as.character(col))
    }
  }
  return(x)
}

#' Will generate a new taxonomic table where the SGB mice names have been replaced with more 
#' interpretable names (prefaced by the taxonomic order they come from). 
#'
#' @param phy_seq_obj A phyloseq object with taxonomic and otu information
#' @export
#' @name clean_SGB_genus
#'
clean_SGB_genus <- function(phy_seq_obj){
  colnames(phyloseq::tax_table(phy_seq_obj)) <-  tolower(phyloseq::rank_names(phy_seq_obj))
  phyloseq::tax_table(phy_seq_obj) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("tax_name_temp") %>%
    dplyr::group_by(phylum) %>%
    dplyr::mutate(
      phylum = paste0("p__", phylum),
      class = ifelse(is.na(class), paste0(phylum, " ", dplyr::row_number()), paste0("c__", class)),
      order = ifelse(is.na(order), paste0(class, " ", dplyr::row_number()), paste0("o__", order)),
      family = ifelse(is.na(family), paste0(order, " ", dplyr::row_number()), paste0("f__", family))) %>% 
    dplyr::mutate(genus = ifelse(startsWith(genus, "GGB"), paste0(family, " ", dplyr::row_number()), ifelse(grepl("^[gfocp]__.+", genus), genus, paste0("g__", genus)))) %>%
    dplyr::mutate(genus = ifelse(startsWith(family, "f__FGB"), paste0(order, " ", dplyr::row_number()), genus)) %>%
    dplyr::mutate(genus = ifelse(startsWith(order, "o__OFG"), paste0(class, " ", dplyr::row_number()), genus)) %>%
    dplyr::mutate(genus = ifelse(startsWith(class, "c__CFG"), paste0(phylum, " ", dplyr::row_number()), genus)) %>%
    dplyr::ungroup() %>%
    tibble::column_to_rownames("tax_name_temp") %>%
    as.matrix() %>% 
    phyloseq::tax_table()
}
#' Will generate a new taxonomic table where the genus names are tidied up
#' interpretable names (prefaced by the taxonomic order thy come from). 
#' This function should be used for cases where the taxonomy is based on 
#' refseq blast hits, such as our amplicon pipeline,
#'
#' @param phy_seq_obj A phyloseq object with taxonomic and otu information
#' @export
#' @name clean_refseq_genus
#'
clean_refseq_genus <- function(phy_seq_obj){
  required_ranks = c("genus", "family", "order", "class", "phylum")
  for (r in required_ranks){
    if (!r %in% tolower(phyloseq::rank_names(phy_seq_obj))){
      stop(paste(r, " must be a tax rank of the phyloseq object.\nCurrent ranks are:", paste(phyloseq::rank_names(phy_seq_obj), sep="", collapse = ", ")))
    }
  }
  phyloseq::tax_table(phy_seq_obj) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("tax_name_temp") %>%
    dplyr::group_by(phylum) %>%
    dplyr::mutate(
      phylum = paste0("p__", phylum),
      class = ifelse(is.na(class), paste0(phylum, " ", dplyr::row_number()), paste0("c__", class)),
      order = ifelse(is.na(order), paste0(class, " ", dplyr::row_number()), paste0("o__", order)),
      family = ifelse(is.na(family), paste0(order, " ", dplyr::row_number()), paste0("f__", family)),
      genus = ifelse(is.na(genus), paste0(family, " ", dplyr::row_number()), ifelse(grepl("^[gfocp]__.+", genus), genus, paste0("g__", genus)))) %>% 
    dplyr::ungroup() %>%
    tibble::column_to_rownames("tax_name_temp") %>%
    as.matrix() %>% 
    phyloseq::tax_table()
}


#' Pad MRNs with zeros to avoid the conflicts between dbs storing mrns as characters vs as integers.
#'
#' @param mrn string or integer of identifier
#' @export
#' @name padMRN
#' @examples
#' mrns <- c(1, "42")
#' padMRN(mrns)
#' #  "00000001" "00000042"
#'
padMRN <- base::Vectorize(USE.NAMES = FALSE, function(mrn) {
  if (is.character(mrn)) {
    tmp <- suppressWarnings(as.integer(mrn))
    if (is.na(tmp)) {
      warning("mrn could not be converted to an integer to pad with zeros")
      return(mrn)
    } else {
      mrn <- as.integer(mrn)
    }
  }
  if (nchar(as.character(mrn)) > 8) warning("mrn integer has more than the maximum expected 8 digits")
  return(stringr::str_pad(as.character(mrn), 8, "left", "0"))
})

#' Make a phyloseq object from some metadata with sampleids
#'
#' @param metadata dataframe with at least one column sampleids; everything else is treated as sample info
#' @param sampleid_col name of the column containing sampleids
#' @param skip_seqs if true, sequences from asv_sequences_ag table will be included in phyloseq object.  Not enabled by default
#' @export
#' @name vdb_make_phylo
#' @examples
#' \dontrun{
#' vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
#' }
vdb_make_phylo <- function(metadata, sampleid_col = "sampleid", skip_seqs = TRUE) {
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame")
  metadata <- process_metadata(metadata, sampleid_col)
  # Make a phyloseq object from some metadata with sampleids
  # TODO add a phylogeny
  assert_db_connected()
  print("getting counts")
  counts <- get_counts_subset(unique(data.frame(metadata)[, sampleid_col]))

  print("creating otu table")
  tab <- counts %>%
    dplyr::select(asv_key, sampleid, count) %>%
    tidyr::pivot_wider(names_from = "asv_key", values_from = "count", values_fill = 0) %>%
    tibble::column_to_rownames("sampleid")
  ot <- phyloseq::otu_table(tab,
    taxa_are_rows = FALSE
  )
  if (!"asv_annotation_blast_ag" %in% ls()) {
    print("getting ASV annotations")
    get_table_from_database("asv_annotation_blast_ag")
  }

  print("creating ASV taxonomy table")
  tax <- asv_annotation_blast_ag[asv_annotation_blast_ag$asv_key %in% unique(counts$asv_key), ] %>%
    dplyr::select(-key, -uploaded_date, -blast_pass) %>%
    dplyr::rename(any_of(c(order = "ordr"))) %>%  # fix typo if exists
    tibble::column_to_rownames("asv_key") %>%
    as.matrix() %>%
    phyloseq::tax_table()

  print("creating sample_data")
  samp <- phyloseq::sample_data(metadata %>% tibble::column_to_rownames(sampleid_col))
  if (!skip_seqs) {
    if (!"asv_sequences_ag" %in% ls()) {
      print("getting ASV sequences")
      get_table_from_database("asv_sequences_ag")
    }

    print("making DNAStringSeq")
    seqs <- asv_sequences_ag[asv_sequences_ag$asv_key %in% unique(counts$asv_key), c("asv_key", "asv_sequence")]
    dss <- Biostrings::DNAStringSet(x = seqs$asv_sequence)
    names(dss) <- seqs$asv_key

    print("constructing phyloseq object")
    return(phyloseq::phyloseq(ot, tax, samp, phyloseq::refseq(dss)))
  } else {
    print("constructing phyloseq object (without sequences)")
  }
  return(phyloseq::phyloseq(ot, tax, samp))
}



get_metaphlan_analyses <- function(con, analysis_ids, schema="public") {
  raw_results <- get_subset_pg_df("mgx_metaphlan", "ia_id", analysis_ids, schema=schema) %>%
    dplyr::filter(grepl("UNCLASSIFIED", clade_name) | grepl(".*\\|s__.*", clade_name)) %>%
    dplyr::filter(!grepl(".*t__.*", clade_name)) %>%
    dplyr::mutate(clade_name = ifelse(clade_name == "UNCLASSIFIED", "k__UNCLASSIFIED", clade_name))
  wide_results <- raw_results %>%
    dplyr::select(ia_id, clade_name, relative_abundance) %>%
    tidyr::pivot_wider(names_from = ia_id, values_from = relative_abundance, values_fill = 0)
  md <- raw_results %>%
    dplyr::filter(clade_name == "k__UNCLASSIFIED") %>%
    dplyr::mutate(estimated_mapped = nreads_input - estimated_number_of_reads_from_the_clade) %>%
    dplyr::select(ia_id, mpa_version, nreads_input, estimated_mapped)
  return(list(wide_results, md))
}

#' Given sample-ids - get the associated experiments and analyses associated with those experiments (optionally filtering by project or app_id)
#'
#' @param sample_ids list of sample ids
#' @param verbose output additional logging information (default = F)
#' @param app_id isabl app identifier (optional)
#' @param proj_id isabl project to filter by
#' @param schema name of schema to query; public for production tables, test_tables for dev, and main when using sqlitedb
#' @param allow_excluded Attempt to include analyses tagged for exclusion
#' @name get_sample_isabl_info
#' @export
#' @examples
#' # assuming test_metadata is a dataframe with a column containing sample IDs called samplid
#' \dontrun{
#' get_sample_isabl_info(list(db[, sample_id_col]), app_id = 66)
#' }
get_sample_isabl_info <- function(sample_ids, verbose=FALSE, app_id=NA, proj_id=NA, schema = "public", allow_excluded=FALSE){
  assert_db_connected()
  db_samples <- get_subset_pg_df("isabl_api_sample", "identifier", sample_ids, schema = schema)
  if (verbose) print(paste("Identified ", nrow(db_samples), "samples of the ", length(sample_ids), " samples requested"))

  db_experiments <- get_subset_pg_df("isabl_api_experiment", "sample_id", list(db_samples$id), schema = schema)
    # If a project id was provided - only return those experiments in the provided project
  if (!is.na(proj_id)){
    db_experiments <- db_experiments %>%
      dplyr::left_join(dplyr::tbl(psql_con, paste(schema, "isabl_api_experiments_projects", sep = ".")) %>%
                       dplyr::select(experiment_id, project_id),
              by=c("id" = "experiment_id")) %>% 
      dplyr::filter(project_id == proj_id) %>%
      dplyr::select(-project_id)
    if (verbose) print(paste("Identified ", nrow(db_experiments), " experiments in Project", proj_id, "associated with those samples"))
  } else{
    if (verbose) print(paste("Identified ", nrow(db_experiments), " experiments associated with those samples"))
  }

  db_analysis_targets <- get_subset_pg_df("isabl_api_analysis_targets", "experiment_id", list(db_experiments$id), schema = schema)
  db_applications <- get_subset_pg_df("isabl_api_application", "id", c(app_id), schema = schema)

  db_analyses <- get_subset_pg_df("isabl_api_analysis", "id", list(db_analysis_targets$analysis_id), schema = schema)
  if(!is.na(app_id)){
    db_analyses <- db_analyses %>%
      dplyr::filter(application_id == app_id) 
  }
  print(db_analyses$exclusion_reason)
  if (!allow_excluded){
    excluded = db_analyses %>% dplyr::filter(exclusion_reason != "" & !is.na(exclusion_reason))
    print(paste0("Dropping ", nrow(excluded), " analyses tagged for exclusion:", paste0(excluded$id, collapse = ",", sep="")))
    db_analyses <- db_analyses %>% dplyr::filter(exclusion_reason == "" | is.na(exclusion_reason))
  }
  if (nrow(db_analyses) == 0){
    if(is.na(app_id)){
      stop(paste("No analyses were found for the samples provided.  These are the samples in isabl:", list(db_samples$id)))
    } else{
      stop(paste("Application", db_applications[db_applications$id == app_id, ][["name"]], "(", db_applications[db_applications$id == app_id, ][["version"]], ") has not been run on any of these samples"))
    }
  }

  db_all_joined <- data.frame(analysis_id = db_analyses$id) %>%
    dplyr::left_join(db_analysis_targets %>% dplyr::select(-id), by = dplyr::join_by(analysis_id)) %>%
    dplyr::left_join(db_experiments %>% dplyr::select(id, identifier, sample_id) %>% dplyr::rename(experiment_identifier = identifier), by = c("experiment_id" = "id")) %>%
    dplyr::left_join(db_samples %>% dplyr::select(id, identifier) %>% dplyr::rename(sample_identifier = identifier), by = c("sample_id" = "id")) 

  return(list(db_all_joined, db_samples, db_experiments, db_analyses, db_analysis_targets))
}

#' Get all analyses for a set of samples, optionally filering by application
#'
#' @param sampleids vector of sample identifiers
#' @param app_id isabl app identifier
#' @param proj_id isabl project to filter by
#' @export
#' @name get_isabl_analyses
#' @examples
#' \dontrun{
#' connect_database()
#' sampleids <- c("2711D", "2711E", "2711F", "2711G", "2711H", "2711I", "2954A", 
#'               "2954B", "2954C", "2954D", "3384D", "3384E", "3384F", "3384H")
#'  get_isabl_analyses(sampleids = sampleids, app_id=43)
#'  }
get_isabl_analyses <- function(sampleids, app_id=NA, proj_id=NA){
  res <- get_sample_isabl_info(sample_ids=sampleids, app_id=app_id, proj_id=proj_id)
  return(res[[4]])
}

#' Make a phyloseq object from some metadata with sampleids
#'
#' @param metadata dataframe with at least one column sampleids; everything else is treated as sample info
#' @param sampleid_col name of the column containing sampleids
#' @param app_id isabl app identifier
#' @param verbose output additional logging information
#' @param choose_max_experiment if multiple analyses are found, select the one utilixing the most experiemnts
#' @param testing use test tables if TRUE
#' @name vdb_make_phylo_mgx
#' @export
#' @examples
#' # assuming test_metadata is a dataframe with a column containing sample IDs called samplid
#' \dontrun{
#' vdb_make_phylo_mgx(test_metadata, sampleid_col = "sampleid")
#' }
vdb_make_phylo_mgx <- function(metadata, sampleid_col = "sampleid", app_id = 66, verbose = FALSE, choose_max_experiment = FALSE, testing = FALSE) {
  
  # If testing we will use the test_schema rather than the real "public" schema.
  schema = "public"
  if (testing){ 
    if (class(psql_con)[1] == "SQLiteConnection"){
      schema = "main"
    } else{ 
      schema = "test_tables"
    }
  }

  # Input checking:
  #------------------------------------------------------------------------------------------------------

  # Currently - this function only supports the Biobakery app:
  db_applications <- get_subset_pg_df("isabl_api_application", "id", c(app_id), schema = schema)
  if (!any(grepl("Biobakery", db_applications["name"]))) {
    stop(paste0("Currently only Biobakery/Metaphlan based make_phylo_mgx is supported. ",
    "Please open a ticket on Github if you need extended functionality! "))
  }

  # check that the sampleid_col is actually in metadata:
  if (!any(sampleid_col %in% colnames(metadata))){
    if (verbose){
      print(paste("Sample id col:", sampleid_col, 
      "not in colnames of metadata supplied:",
      colnames(metadata),"check for typos and try again!"))
    }
    stop("sampleid_col not in colnames of metadata supplied.")
  }

  # Pulling together tables from relevant sources section:
  #------------------------------------------------------------------------------------------------------

  # process metadata, and use it to fetch isabl tables related to these sample ids. 
  metadata <- process_metadata(metadata, sampleid_col)
  isabl_res <- get_sample_isabl_info(
    list(metadata[, sampleid_col]), 
    verbose=verbose, app_id=app_id, 
    schema = schema)

  # pull out the joined db of all samples, experiments and analyses corresponding to these sample ids. 
  db_all_joined <- isabl_res[[1]] %>%
    dplyr::group_by(analysis_id) %>%
    dplyr::mutate(num_exps_per_analysis =  dplyr::n()) %>%
    dplyr::mutate(sample_level = dplyr::n() > 1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(num_experiments_per_sample =  length(unique(experiment_id))) %>%
    dplyr::mutate(num_analyses_per_sample =  length(unique(analysis_id))) %>%
    dplyr::ungroup()

    # fetch the metaphlan results for these analyses and parse out the results and metadata:
  res <- get_metaphlan_analyses(con = psql_con, analysis_ids = unique(db_all_joined$analysis_id), schema = schema)
  results <- res[[1]]
  md <- res[[2]]

  # Tests/sanity checks on the analyses/results section:
  #------------------------------------------------------------------------------------------------------
  
  # Check and (optionally fix) the samples where multiple passes of the analysis have been run on different 
  # experiments associated with the sample. (This can happen for example if the analysis was run all all 
  # sequences of a sample, but then the sample was later resequenced, and a new copy of the analysis was run 
  # including the newly sequenced data.) By setting choose_max_experiment = T this code will choose (one of) 
  # the analyses with the highest number of other experiments in that analysis.  If it is not set to true, this
  # code will error and request that the error be fixed in isabl directly. 
  if (any(db_all_joined$num_analyses_per_sample > 1)){
    if(choose_max_experiment){
      if (verbose) print(paste0("Since choose_max_experiment is true,",
        " will select the analysis which has the most experiments associated",
        " with it to proceed."))
        db_all_joined <- db_all_joined %>%
          dplyr::group_by(sample_identifier) %>%
          dplyr::mutate(max_analysis = analysis_id[which.max(num_exps_per_analysis)]) %>% # which.max just chooses one on a tie
          dplyr::ungroup() %>%
          dplyr::filter(analysis_id == max_analysis)
        md <- md %>%
          dplyr::filter(ia_id %in% db_all_joined$max_analysis)

    }
    else{
      print(db_all_joined %>% dplyr::filter(num_exps_per_analysis > 1))
      stop(paste0("The same experiments are in multiple analyses. ", 
        "This can happen when someone runs an app both per-sample and per-experiment. ",
        "Resolve this in Isabl to avoid confusion, or rerun this code with ",
        "choose_max_experiment = TRUE.  The problem samples are above."
        ))
    }
  }

  md <- md %>%
    dplyr::rename(analysis_id = ia_id) %>%
    dplyr::full_join(db_all_joined, by = dplyr::join_by(analysis_id))
  # check if there a nasty mix of sample-level and experiment-level analyses. We don't want these.
  # previous version of isabl_microbiome_apps permitted experiment-level runs.
  problematic_analyses <- md %>% dplyr::filter(!sample_level, num_experiments_per_sample > 1) 
  if(nrow(problematic_analyses) > 0){
    print("Problematic analyses found!")
    print(problematic_analyses)
    warning(paste0("These analyses are experiment-level but should be deleted and re-run at the sample-level. They will be removed for now, but need to be fixed at the source"))
    md <- md %>% dplyr::filter(!analysis_id %in% problematic_analyses$analysis_id)
  }

  # Create phyloseq section:
  #------------------------------------------------------------------------------------------------------
  if (verbose) print("creating otu table")
  ot <- phyloseq::otu_table(results %>% tibble::column_to_rownames("clade_name"), taxa_are_rows = TRUE)
  tax_tab <- dplyr::bind_rows(lapply(gsub("\\|", ";", rownames(ot)), phyloseq::parse_taxonomy_qiime)) %>% as.matrix()
  rownames(tax_tab) <- rownames(ot)

  if (verbose) print("creating sample_data")
  abbreviated_md <- md %>%
    dplyr::select(-c(experiment_id, experiment_identifier)) %>%
    unique()
  
  samp <- metadata %>% dplyr::rename(sample_identifier = dplyr::all_of(sampleid_col)) %>%
    dplyr::full_join(abbreviated_md, by = dplyr::join_by(sample_identifier)) %>%
    tibble::column_to_rownames("sample_identifier") 
  if (any(is.na(samp$analysis_id))){
    warning("Removing the following samples from metadata lacking analyses")
    cat(paste0("\n\n-----------------------\nThese are the samples lacking analyses: \n\n",
    sum(is.na(samp$analysis_id)), " NA analyses removed."))
    print(samp[is.na(samp$analysis_id), ])
    samp <- samp %>% dplyr::filter(!is.na(analysis_id))
  }
  
  sample_names_df  <- db_all_joined %>%
    dplyr::mutate(analysis_id = as.character(analysis_id)) %>%
    dplyr::select(c(analysis_id, sample_identifier)) %>%
    unique() %>%
    dplyr::right_join(data.frame(analysis_ids = phyloseq::sample_names(ot)), by=c("analysis_id" = "analysis_ids"))
  
  phyloseq::sample_names(ot) = data.frame(analysis_id=phyloseq::sample_names(ot)) %>% dplyr::left_join(sample_names_df, by="analysis_id") %>% dplyr::pull(sample_identifier)

  if (verbose) print("contructing phyloseq object")
  return(phyloseq::phyloseq(ot, phyloseq::tax_table(tax_tab), phyloseq::sample_data(samp)))
}


#' Helper function that speeds up interactions with the postgres database by only
#' fetching those parts of the table which contain a list of desired IDs. 
#' @param table_name the name of the table to do the subsetting on.
#' @param table_id_name the name of the id in the table to perform filtering with.  
#' @param ids list of strings: ids to use in the query
#' @param con - the postgres connection object, defaults to psql_con
#' @param schema name of schema to query; public for production tables, test_tables for dev

#' @export
#' @name get_subset_pg_df
get_subset_pg_df <- function(table_name, table_id_name, ids, con = NA, schema = "public") {
  if(is.na(con)){
    con = psql_con
  }
  ids_str = paste("'", unlist(ids), "'", sep = "", collapse = ", ")
  query_str = paste0("SELECT * FROM ", schema, ".", table_name, " WHERE ", table_id_name, " in (", ids_str, ")")

  
  print(query_str)
  return(RPostgres::dbGetQuery(con, query_str))
}

#' For a given project number, get all the samples assocoiated with it - optionally
#' filtering by those samples which have successful analyses run through a particular 
#' application number. 
#'
#' @param project the project number to fetch analysis ids for
#' @param app_id the number of the application to check samples have a successful analysis for. 
#' @param verbose if messages should be printed about what is being run/found. 
#' @export
#' @name get_project_samples
#' @examples
#' # ie to get all the Human Biobakery results for project 32: 
#' \dontrun{
#' ia_ids = get_project_samples(project = 32, app_id = 66, verbose = T)
#' }
get_project_samples <- function(project, app_id = NA, verbose = FALSE) {
  assert_db_connected()
  exp_query = paste0("SELECT * FROM isabl_api_experiments_projects WHERE project_id in (", project, ")")
  exp_ids <- RPostgres::dbGetQuery(psql_con, exp_query)
  if (verbose) print(paste("Found ", nrow(exp_ids), " experiments associated with project ", project))

  # If passed an application ID we will only return those samples which have an analysis of that type performed. 
  # Note - it doesn't filter on whether or not said analysis was successful though!
  if(!is.na(app_id)){
    exp_ids <- exp_ids %>%
      dplyr::right_join(get_subset_pg_df("isabl_api_analysis_targets", "experiment_id", list(exp_ids$experiment_id)), by = "experiment_id")
    if (verbose) print(paste("Found ", nrow(exp_ids), " analyses associated with these experiments."))
    exp_ids <- exp_ids %>%
      dplyr::inner_join(get_subset_pg_df("isabl_api_analysis", "id", list(exp_ids$analysis_id)), by = c("analysis_id" = "id")) %>% 
      dplyr::filter(application_id == app_id)
    if (verbose) print(paste("After filtering, found  ", nrow(exp_ids), " experiments which are associated with an experiment from application: ", app_id))
  }
  
  exp_ids <- exp_ids %>%
    dplyr::left_join(get_subset_pg_df("isabl_api_experiment", "id", list(exp_ids$experiment_id)), by = c("experiment_id" = "id")) %>%
    dplyr::rename(experiment_identifier = identifier)
  sample_ids <- get_subset_pg_df("isabl_api_sample", "id", list(exp_ids$sample_id)) %>%
    dplyr::rename(sample_identifier = identifier) %>%
    dplyr::rename(samples_table_id = id)
  exp_ids <- exp_ids %>%
    dplyr::left_join(sample_ids, by = c("sample_id" = "samples_table_id"))
  return(unique(exp_ids$sample_identifier))
}
