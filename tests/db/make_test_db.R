library(vdbR)
library(tidyverse)
get_table_from_database("frozen_set_ag")
get_table_from_database("asv_annotation_blast_color_ag")
get_table_from_database("asv_annotation_blast_ag")
set.seed(123)
tmp <- frozen_set_ag %>% dplyr::filter(institution=="MSK_allo") %>% dplyr::sample_n(size = 50) %>% 
  mutate(institution = sample(c("A", "B"), size=n(), replace=T))  %>% 
  select(-hct, -patient_id, -uploaded_date, -key)
write.csv(tmp,"tests/tmp.samples.csv", row.names = FALSE)

counts <- get_counts_subset(tmp$sampleid)
counts %>% select(-key, -uploaded_date, -count_relative)   %>% write.csv("tests/tmp.counts.csv",  row.names = FALSE)
asv_annotation_blast_ag %>% filter(asv_key %in% counts$asv_key) %>%  write.csv("tests/tmp.anno.csv",  row.names = FALSE)

asv_annotation_blast_color_ag %>% filter(asv_key %in%counts$asv_key) %>% dplyr::select(-key, -uploaded_date) %>% write.csv("tests/tmp.colors.csv",  row.names = FALSE)


if(file.exists( "inst/extdata/db.sqlite")) system("rm  inst/extdata/db.sqlite")
system("sqlite3 inst/extdata/db.sqlite < tests/db/populate_test_db.sql ")

# test with
# psql_con = DBI::dbConnect(RSQLite::SQLite(), "tests/db.sqlite")
# tbl <- get_counts_subset(tmp$sampleid)
# or
# vdb_make_phylo(tmp, psql_con = psql_con)
