# vdbR
Set of functions to help accessing microbiota database and and analyzing data

## Install
Install using this command

```
remotes::install_github("vdblab/vdbR")

```


## Useful functions
### Connecting to the database
The database config file is a comma-delimited file as follows:

```
user,pass,host,dbname
me,test123,my.postgres.url.org,my_db_name

```

Your database config file should be placed in one of the following approved locations:

  - `~/dbConfig.txt`
  - `~/.dbConfig.txt`
  - `~/.config/dbConfig.txt`
  - `~/.config/vdbR/dbConfig.txt`


This prevents having to hardcode the path to the config file in your scripts. With your config in one of the approved locations, this connects you to the microbiome database server:
```r
connect_database()
```


### Database Usage

List all the available tables:
```r
list_table_from_database(pattern="*")
```

Get a table from the database (loads it to your global env):
```r
get_table_from_database("qpcr_16s_ag")
dim(qpcr_16s_ag)
# [1] 4748    7
```


## Make a phyloseq ID based on some metadata table
```
test_metadata <- structure(
  list(
    sampleid = c("2133C", "1773A", "1773I", "90.tp.51", 
                 "223D", "2065C", "1773O", "1773M", "2133D", "2065D", "90.tp.97", 
                 "1773E", "2065B", "2065E", "90.tp.57", "2133G", "2065F", "1773G", 
                 "FMT.0113I", "2065G"), 
    age = c("60 and over ", 
            "60 and over ", "60 and over ", "60 and over ", "60 and over ", 
            "40-59", "60 and over ", "60 and over ", "60 and over ", "40-59", 
            "40-59", "60 and over ", "40-59", "40-59", "40-59", "60 and over ", 
            "40-59", "60 and over ", "40-59", "40-59"),
    fakemrn = c("00000001", 
                "00000002", "00000003", "00000004", "00000005", "00000001", "00000002", 
                "00000003", "00000004", "00000005", "00000001", "00000002", "00000003", 
                "00000004", "00000005", "00000001", "00000002", "00000003", "00000004", 
                "00000005")), 
  row.names = c(NA, -20L), class = "data.frame")

vdb_make_phylo(test_metadata, sampleid_col = "sampleid")
```


## Get a dataframe of all Isabl analyses for a list of samples

```
sampleids <- c("2711D", "2711E", "2711F", "2711G", "2711H", "2711I", "2954A")
get_isabl_analyses(sampleids = sampleids, app_id=43)

```

## Get QC metric for Shotgun samples
```
sampleids <-  c("2711D", "2711E", "2711F", "2711G", "2711H", "2711I", "2954A")
preproces_analyses <- get_isabl_analyses(sampleids = sampleids, app_id=43)
# note there are more rows than samples: one per flowcell

dplyr::tbl(psql_con, "mgx_qc") %>% 
   dplyr::filter(ia_id %in% local(preproces_analyses$ia_id)) %>% 
   dplyr::collect()


```

## Make Shotgun Phyloseq Object

```
metadata <-  data.frame(identifier = sampleids)
phy <- vdb_make_phylo_mgx(metadata = metadata, sampleid_col = "identifier", app_id = 66, verbose = FALSE)
phy
phyloseq::plot_bar(phy, fill="Order")

```

We have currently implemented an easy way to get a phyloseq object from metaphlan results stored in the database.  to do this, you need (a) some metadata about what samples you want, and (b) and Isabl application ID.  Each app version gets its own ID, so specifying an application ID ensures that the results will be compatible with eachother. 

For a list of apps, see 
```
connect_database()
dplyr::tbl(psql_con, "isabl_api_application") %>% dplyr::collect()
```



## Scripts

## Development
### Update docs and NAMESPACE
```
devtools::document()
```
### Increment version using:
```
usethis::use_version()
```

### Style package
```
styler::style_pkg(".")
```

### Style tutorials
```
styler::style_dir("data_analysis_tutorial/", filetype = "rmd")
```
### Auto-Name tutorial chunks
```
namer::name_dir_chunks("data_analysis_tutorial/")
```

### Render Tutorials
```
bookdown::render_book("data_analysis_tutorial/")
```

## Admin Tasks
These tasks aren't run by the average user, so they are not exported when loading the package.  Access these functions using the `:::` operator.
### Add a new user
```
library(vdbR)
connect_database()
vdbR:::create_new_postgres_user(user_name = "nwtest", temp_pass = "test123456")
```

### Check to see if user changed their temporary password
```
vdbR:::block_user_who_did_not_change_password(user_name = "nwtest", temp_pass = "test123456")
```

### Updating database tables
TBD


# History
Mirrored from Antonio Gomes's original repo at github.mskcc.org/gomesa/vdbr, and added various scripts sporadically developed over the years for database admin, processing data, and more.
Moved from vdblabinternal to the Microbiome organization
