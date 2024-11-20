# Apr/20/2020
#
# Making package of connecting to database
# Install: remotes::install_git("https://github.mskcc.org/gomesa/ag.gutbug.git")
#
# Aug/09/2019
# We had the `maximum connections` capacity reached and couldn't connect anymore. Many of the connections of from user `postgres`.
# I solved it by forcing to terminate all connection:
# SELECT *, pg_terminate_backend(pid) FROM pg_stat_activity;
#
# After connection, use this to see list of tables:
#     dbListTables(con)
# Use this to load a table
# My_table = dbReadTable(psql_con,TABLE_NAME)
#
#
# Use the this function in the attached script to connect:
#   connect_database(username, password);
# At your first usage, change the password:
#   dbSendQuery(conn = psql_con,"alter role YOUR_USERNAME password 'YOUR_PASSWORD'”);
# Use this function to list table from database:
#   list_table_from_database()
# You can use this function to load a table:
#   get_table_from_database()
#






approved_locations <- file.path(Sys.getenv("HOME"), c(
  "dbConfig.txt",
  ".dbConfig.txt",
  ".config/dbConfig.txt",
  ".config/vdbR/dbConfig.txt"
))

check_config_location <- function(config_file) {
  if (missing(config_file)) {
    if (any(file.exists(approved_locations))) {
      config_file <- approved_locations[file.exists(approved_locations)][[1]]
    } else {
      stop("No config file found in any of the approved locations:\n", paste(approved_locations, collapse = "\n"))
    }
  } else {
    if (!file.exists(config_file)) stop(paste("Config file ", config_file, " does not exist"))
  }
  if (!config_file %in% approved_locations) {
    warning(paste0("Please move your config file to one of the following locations to ensure portability:\n  ", paste0(approved_locations, collapse = "\n  ")))
  }
  return(config_file)
}




#' connect_database
#' @name connect_database
#' @import data.table
#' @param config_file file with login info to connect to database
#' @param bundled conected to test sqlite database
#' @export


connect_database <- function(config_file, bundled=FALSE) {
  if (bundled){
    psql_con = DBI::dbConnect(
      RSQLite::SQLite(), 
      system.file("extdata", "db.sqlite", package = "vdbR") 
    )
    assign("psql_con",
           psql_con,
           envir = .GlobalEnv
    )
    return()
  }
  config_file <- check_config_location(config_file = config_file)

  config <- data.table::fread(config_file, nrows = 1)

  drv <- RPostgres::Postgres()

  # loads the PostgreSQL driver
  # creates a connection to the postgres database
  # note that "con" will be used later in each connection to the database
  psql_con <<- RPostgres::dbConnect(drv,
    dbname = config$dbname,
    host = config$host, port = 5432,
    user = config$user, password = config$pass
  )

  if (config$pass == "test123456") {
    print("You are using the default password `test123456`")
    print("Please, update your password: DO NOT USE YOUR MSK PASSWORD!")
    new_pass <- readline(prompt = "type your new password:")
    pw_update_line <- sprintf(
      "alter role %s password '%s'",
      config$user,
      new_pass
    )
    RPostgres::dbSendQuery(conn = psql_con, pw_update_line)
    config$pass <- new_pass
    write.table(config, config_file, row.names = F, quote = F, sep = ",")
    print("Your password was updated! Config file was changed accordingly")
  }

}

#' assert_db_connected
#' @name assert_db_connected
#' @param con string of global variable containing the connection, defaults to psql_con created by connect_database
#' @export

assert_db_connected <- function(con="psql_con") {
  # TODO: this relies on scoping, which complicates testing unless you set global variables in tests
  #  if (!con %in% ls(envir =parent.frame())){
  if (!exists(con)){
      # perhaps we should include the test db as package data instead
    stop(paste0(con, " not found; please connect to production database by running `connect_database()` or to the test database with `connect_database(bundled=TRUE)`"))
  }
    
}
 
