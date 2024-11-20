# 

#' Querying database
#' @name get_table_from_database
#' @export
#' @param table_name name of table to be loaded.
#' @param pre_filter If True for asv_counts_ag and asv_alpha_diversity_ag,it selects only run of highest coverage in case of multiple runs for same sample.
#' @param pre_filter If True for samples_castori_ag, it filters out sampleids that appear in samples_castori_blacklist.
#' @param pattern only tables that match \emph{pattern} will be returned.
#' @param sampleid_subset Only sampleids in this vector will be returned.
#'
get_table_from_database <- function(table_name, pre_filter = TRUE) {
  assert_db_connected()
  if (class(psql_con)[1] == "SQLiteConnection"){
    tb <- RSQLite::dbReadTable(psql_con, table_name)
  } else{
    tb <- RPostgres::dbReadTable(psql_con, table_name)
  }
  tb <- data.table::data.table(tb)
  assign(table_name, tb, envir = .GlobalEnv)

  pre_filter_set <- c("asv_counts_ag", "asv_alpha_diversity_ag", "samples_castori_ag")
  if (pre_filter & table_name %in% pre_filter_set) {
    if (table_name == "asv_counts_ag") {
      asv_counts_ag <- tb
      oligos_id_filtered <- unique(asv_counts_ag[count_total > 200][order(sampleid, -count_total)][!duplicated(sampleid)]$oligos_id)
      asv_counts_ag <- asv_counts_ag[count_total > 200][order(sampleid, -count_total)][oligos_id %in% oligos_id_filtered]
      asv_counts_ag[, count_relative := count / count_total]
      asv_counts_ag$filtered_for_highest_coverage_run <- TRUE
      asv_counts_ag <<- asv_counts_ag
      print(sprintf(
        "table %s is loaded and filtered for duplicates. Only the replicate of highest coverage is retained.",
        table_name
      ))
    } else if (table_name == "asv_alpha_diversity_ag") {
      asv_alpha_diversity_ag <- tb
      oligos_id_filtered <- unique(asv_alpha_diversity_ag[count_total > 200][order(sampleid, -count_total)][!duplicated(sampleid)]$oligos_id)
      asv_alpha_diversity_ag <- asv_alpha_diversity_ag[count_total > 200][order(sampleid, -count_total)][oligos_id %in% oligos_id_filtered]
      asv_alpha_diversity_ag$filtered_for_highest_coverage_run <- TRUE
      asv_alpha_diversity_ag <<- asv_alpha_diversity_ag
      print(sprintf(
        "table %s is loaded and filtered for duplicates. Only the replicate of highest coverage is retained.",
        table_name
      ))
    } else if (table_name == "samples_castori_ag") {
      samples_castori_ag <- tb
      samples_castori_blacklist <- RPostgres::dbReadTable(psql_con, "samples_castori_blacklist")
      samples_castori_ag <- samples_castori_ag[!samples_castori_ag$sampleid %in% samples_castori_blacklist$sampleid, ]
      samples_castori_ag <<- samples_castori_ag
      print(sprintf(
        "table %s is loaded and filtered for blacklisted samples.",
        table_name
      ))
    }
  }
}

#' @rdname get_table_from_database
#' @export
list_table_from_database <- function(pattern = NULL) {
  tb <- RPostgres::dbListTables(psql_con)

  if (!is.null(pattern)) {
    tb <- tb[grep(pattern, tb, ignore.case = T)]
  }
  return(sort(tb))
}

#' @rdname get_table_from_database
#' @export
get_counts_subset <- function(sampleid_subset = NULL, pre_filter = T) {
  if (is.null(sampleid_subset)) {
    return(NULL)
  } else {
    sampleid_subset_str <- paste(sprintf("'%s'", sampleid_subset), collapse = ",")
    where_query <- sprintf(" where  sampleid in (%s)", sampleid_subset_str)

    query <- sprintf(
      "select * from asv_counts_ag %s",
      where_query
    )
    # check if using the sqlite test db
    if (class(psql_con)[1] == "SQLiteConnection"){
      asv_counts_subset <- RSQLite::dbGetQuery(psql_con, query)
    } else{
      asv_counts_subset <- RPostgres::dbGetQuery(psql_con, query)
    }
    asv_counts_subset <- data.table::as.data.table(asv_counts_subset)

    if (pre_filter) {
      # print(class(asv_counts_subset));
      oligos_id_filtered <- unique(asv_counts_subset[count_total > 200][order(sampleid, -count_total)][!duplicated(sampleid)]$oligos_id)
      asv_counts_subset <- asv_counts_subset[count_total > 200][order(sampleid, -count_total)][oligos_id %in% oligos_id_filtered]
      asv_counts_subset[, count_relative := count / count_total]
    }
    # dim(asv_counts_subset);
    return(asv_counts_subset)
  }
}
