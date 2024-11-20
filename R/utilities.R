#' Ensure that metadata is a dataframe without NA sampleid or duplicated sampleids.
#' cast metadata to dataframe so we can do convenient indexing with metadata[1, samplid_col]
#'
#' @param metadata metadata dataframe
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
#' @param phy_seq_obj
#' @export
#' @name clean_SGB_genus
#'
clean_SGB_genus <- function(phy_seq_obj){
  colnames(phyloseq::tax_table(phy_seq_obj)) <-  tolower(phyloseq::rank_names(phy_seq_obj))
  phyloseq::tax_table(phy_seq_obj) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("tax_name_temp") %>%
    dplyr::group_by(phylum) %>%
    dplyr::mutate(genus = ifelse(startsWith(genus, "GGB"), paste0("f__", family, " ", dplyr::row_number()), ifelse(grepl("^[gfocp]__.+", genus), genus, paste0("g__", genus)))) %>%
    dplyr::mutate(genus = ifelse(startsWith(family, "FGB"), paste0("o__", order, " ", dplyr::row_number()), genus)) %>%
    dplyr::mutate(genus = ifelse(startsWith(order, "OFG"), paste0("c__", class, " ", dplyr::row_number()), genus)) %>%
    dplyr::mutate(genus = ifelse(startsWith(class, "CFG"), paste0("p__", phylum, " ", dplyr::row_number()), genus)) %>%
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
#' @param phy_seq_obj
#' @export
#' @name clean_refseq_genus
#'
clean_refseq_genus <- function(phy_seq_obj){
  newranks <-  tolower(phyloseq::rank_names(phy_seq_obj))
  # for historical purposes:
  if ("ordr" %in% newranks){
    newranks[which(newranks == "ordr")] <- "order"
  } 
  colnames(phyloseq::tax_table(phy_seq_obj)) <-  newranks

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
      order = ifelse(is.na(order), paste0("c__", class, " ", dplyr::row_number()), order),
      family = ifelse(is.na(family), paste0("o__", order, " ", dplyr::row_number()), family),
      genus = ifelse(is.na(genus), paste0("f__", family, " ", dplyr::row_number()), ifelse(grepl("^[gfocp]__.+", genus), genus, paste0("g__", genus)))) %>% 
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
#' @example
#' mrns <- c(1, "42")
#' padMRN(mrns)
#' [1] "00000001" "00000042"
#'
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
#' @example
#' # assuming test_metadata is a dataframe with a column containing sample IDs called samplid
#' vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
#'
#'

vdb_make_phylo <- function(metadata, sampleid_col = "sampleid", skip_seqs = TRUE, psql_con=NULL) {
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame")
  metadata <- process_metadata(metadata, sampleid_col)
  #' Make a phyloseq object from some metadata with sampleids
  # TODO add a phylogeny
  print("connecting to database")
  if (missing(psql_con)) connect_database()

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

    print("contructing phyloseq object")
    return(phyloseq::phyloseq(ot, tax, samp, phyloseq::refseq(dss)))
  } else {
    print("contructing phyloseq object (without sequences)")
  }
  return(phyloseq::phyloseq(ot, tax, samp))
}



get_metaphlan_analyses <- function(con, analysis_ids) {
  raw_results <- dplyr::tbl(psql_con, "mgx_metaphlan") %>%
    dplyr::filter(ia_id %in% analysis_ids) %>%
    dplyr::collect() %>%
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


#' Make a phyloseq object from some metadata with sampleids
#'
#' @param metadata dataframe with at least one column sampleids; everything else is treated as sample info
#' @param sampleid_col name of the column containing sampleids
#' @param app_id isabl app identifier
#' @param verbose output additional logging information
#' @export
#' @name vdb_make_phylo_mgx
#' @example
#' # assuming test_metadata is a dataframe with a column containing sample IDs called samplid
#' vdb_make_phylo_mgx(test_metadata, sampleid_col = "sampleid")
vdb_make_phylo_mgx <- function(metadata, sampleid_col = "sampleid", app_id = 66, verbose = FALSE) {
  
  metadata <- process_metadata(metadata, sampleid_col)
  if (verbose) print("connecting to database")
  connect_database()
  db_samples <- dplyr::tbl(psql_con, "isabl_api_sample") %>%
    dplyr::collect() %>%
    dplyr::filter(identifier %in% metadata[, sampleid_col])
  if (verbose) print(paste("Identified ", nrow(db_samples), "samples of the ", nrow(metadata), " samples requested"))

  db_experiments <- dplyr::tbl(psql_con, "isabl_api_experiment") %>%
    dplyr::filter(sample_id %in% local(db_samples$id)) %>%
    dplyr::collect()
  if (verbose) print(paste("Identified ", nrow(db_experiments), " experiments associated with those samples"))

  db_applications <- dplyr::tbl(psql_con, "isabl_api_application") %>% dplyr::collect()
  db_analysis_targets <- dplyr::tbl(psql_con, "isabl_api_analysis_targets") %>%
    dplyr::filter(experiment_id %in% local(db_experiments$id)) %>%
    dplyr::collect()

  db_analyses <- dplyr::tbl(psql_con, "isabl_api_analysis") %>%
    dplyr::filter(application_id == app_id) %>%
    dplyr::filter(id %in% local(db_analysis_targets$analysis_id)) %>%
    dplyr::collect()
  if (nrow(db_analyses) == 0) stop(paste("Application", db_applications[db_applications$id == app_id, ][["name"]], "(", db_applications[db_applications$id == app_id, ][["version"]], ") has not been run on any of these samples"))
  targets_per_analysis <- db_analysis_targets %>%
    dplyr::filter(analysis_id %in% db_analyses$id) %>%
    dplyr::group_by(analysis_id) %>%
    dplyr::summarize(n_targets = dplyr::n_distinct(experiment_id)) %>%
    dplyr::pull(n_targets)


  if (!all(targets_per_analysis == 1)) {
    if (dplyr::n_distinct(targets_per_analysis) > 1) {
      if (verbose) print(paste("Analyses are a mix of one or more experiments"))
    } else {
      if (verbose) print(paste("Analyses contain more than one experiments"))
    }
  }


  if (grepl("Biobakery", db_applications[db_applications$id == app_id, ]["name"])) {
    res <- get_metaphlan_analyses(con = psql_con, analysis_ids = db_analyses$id)
    results <- res[[1]]
    md <- res[[2]]

    analysis_ids <- data.frame(analysis_id = as.numeric(colnames(results)[2:ncol(results)]))
    samples_per_analysis_id <- analysis_ids %>%
      dplyr::left_join(db_analysis_targets %>% dplyr::select(-id), by = dplyr::join_by(analysis_id)) %>%
      dplyr::left_join(db_experiments %>% dplyr::select(id, identifier, sample_id) %>% dplyr::rename(experiment_identifier = identifier), by = c("experiment_id" = "id")) %>%
      dplyr::left_join(db_samples %>% dplyr::select(id, identifier), by = c("sample_id" = "id")) %>%
      dplyr::select(-sample_id)
    persample=FALSE
    if (nrow(samples_per_analysis_id) != dplyr::n_distinct(samples_per_analysis_id$identifier)) {
      if (verbose) print("at least one analysis contained more than one experiment  -- returning a per-sample label for those samples")
      persample=TRUE
    }
    if (nrow(samples_per_analysis_id) != dplyr::n_distinct(samples_per_analysis_id$experiment_id)) {
      if (verbose) print(samples_per_analysis_id)
      stop("the same experiments are in multiple analyses. This can happen when someone runs an app both per-sample and per-experiment. Resolve this in Isabl to avoid confusion")
    }
    analysis_ids_tidy <- samples_per_analysis_id %>%
      dplyr::group_by(analysis_id) %>% 
      dplyr::summarize(label=ifelse(dplyr::n_distinct(experiment_identifier) > 1, identifier, experiment_identifier)) %>% 
      tibble::column_to_rownames("analysis_id")
    old_col_names <- colnames(results)[2:ncol(results)]
    colnames(results) <- c("clade_label", analysis_ids_tidy[old_col_names, ])
    md <- md %>%
      dplyr::rename(analysis_id = ia_id) %>%
      dplyr::full_join(samples_per_analysis_id, by = dplyr::join_by(analysis_id)) %>% 
      dplyr::select(-experiment_id) %>%
      dplyr::group_by(dplyr::across(c(-experiment_identifier))) %>% 
      dplyr::summarise(experiments=paste0(experiment_identifier, collapse=",")) %>% 
      dplyr::group_by(analysis_id) %>%
      dplyr::mutate(sample_level = dplyr::n() > 1) %>%
      dplyr::ungroup() %>% 
      dplyr::group_by(identifier) %>% 
      dplyr::mutate(n_experiments = dplyr::n_distinct(experiments)) %>%
      dplyr::ungroup()
  }
  # check if there a nasty mix of sample-level and experiment-level analyses. We don't want these.
  # previous version of isabl_microbiome_apps permitted experiment-level runs.
  problematic_analyses <- md %>% dplyr::filter(!sample_level, n_experiments > 1) 
    
  if(nrow(problematic_analyses) > 0){
    print(problematic_analyses)
    warning(paste0("These analyses are experiment-level but should be deleted and re-run at the sample-level. They will be removed for now, but need to be fixed at the source"))
    md <- md %>% dplyr::filter(!analysis_id %in% problematic_analyses$analysis_id)
  }
  if (verbose) print("creating otu table")
  ot <- phyloseq::otu_table(results %>% tibble::column_to_rownames("clade_label"), taxa_are_rows = TRUE)

  tax_tab <- dplyr::bind_rows(lapply(gsub("\\|", ";", rownames(ot)), phyloseq::parse_taxonomy_qiime)) %>% as.matrix()
  rownames(tax_tab) <- rownames(ot)

  if (verbose) print("creating sample_data")
  samp <- metadata %>% dplyr::rename(identifier = all_of(sampleid_col)) %>%
    dplyr::full_join(md, by = dplyr::join_by(identifier)) %>%
#    group_by(identifier) %>%
#    mutate(identifier = ifelse(persample))
    tibble::column_to_rownames("identifier") 
  
  if (any(is.na(samp$analysis_id))){
    warning("Removing the following samples from metadata lacking analyses")
    print(samp[is.na(samp$analysis_id), ])
    samp <- samp %>% dplyr::filter(!is.na(analysis_id))
  }
  ## we should probably extract the sample identifiers more robustly
  sample_identifiers <- gsub(".*?_(.*)_.*", "\\1", phyloseq::sample_names(ot))
  if(persample & length(sample_identifiers)==length(unique(sample_identifiers))) {
    phyloseq::sample_names(ot) = sample_identifiers
  }else{
    rownames(samp) <-  samp$experiments
  }
  
  if (verbose) print("contructing phyloseq object")
  return(phyloseq::phyloseq(ot, phyloseq::tax_table(tax_tab), phyloseq::sample_data(samp)))
}

#' Get all analyses for a set of samples, optionally filering by application
#'
#' @param sampleids vector of sample identifiers
#' @param app_id isabl app identifier
#' @export
#' @name get_isabl_analyses
#' @example
#'  sampleids <- c("2711D", "2711E", "2711F", "2711G", "2711H", "2711I", "2954A", 
#'               "2954B", "2954C", "2954D", "3384D", "3384E", "3384F", "3384H")
#'  get_isabl_analyses(sampleids = sampleids, app_id=43)
get_isabl_analyses <- function(sampleids, app_id=NA){
  connect_database()
  db_samples <- dplyr::tbl(psql_con, "isabl_api_sample") %>%
    dplyr::select(identifier, id) %>% 
    dplyr::filter(identifier %in% local(sampleids)) %>% 
    dplyr::left_join(dplyr::tbl(psql_con, "isabl_api_experiment") %>%
                       dplyr::select(id, system_id, sample_id, identifier) %>%
                       dplyr::rename(exp_identifier = identifier, experiment_id = id, exp_system_id = system_id),
              by=c("id"="sample_id")) %>% 
    dplyr::rename("sampleid"="id")
  if (is.na(app_id)){
    apps <- dplyr::tbl(psql_con, "isabl_api_application") 
  } else{
    apps <- dplyr::tbl(psql_con, "isabl_api_application")  %>% 
      dplyr::filter(id == app_id) 
  }
  db_analyses <-  apps %>% 
    dplyr::select(id, name, version) %>%
    dplyr::rename(application_id = id ) %>%
    dplyr::left_join(dplyr::tbl(psql_con, "isabl_api_analysis") %>% dplyr::select( id, storage_url, application_id), by ="application_id") %>%
    dplyr::left_join(dplyr::tbl(psql_con, "isabl_api_analysis_targets") %>% dplyr::select(-id), by=c("id"="analysis_id")) %>% 
    dplyr::inner_join(db_samples %>% dplyr::select(identifier, exp_identifier, experiment_id ), by ="experiment_id") %>% 
    dplyr::rename(ia_id = id) 
  return(db_analyses %>% dplyr::collect())
}
#sampleids <- c("2711D", "2711E", "2711F", "2711G", "2711H", "2711I", "2954A", 
#               "2954B", "2954C", "2954D", "3384D", "3384E", "3384F", "3384H")
#get_isabl_analyses(sampleids = sampleids)
