test_that("non-existant config location fails", {
  expect_error(check_config_location(config_file = "/home/turnip"))
})


test_that("weird config name raises warning", {
  expect_warning(check_config_location(config_file = "test-config-and-dbconnect.R"))
})

test_that("non-existant warning", {
  config_file <- approved_locations[file.exists(approved_locations)][[1]]
  expect_silent(check_config_location(config_file = config_file))
})


test_that("missing db connection raises error", {
  expect_error(assert_db_connected(con="turnip"))
})

test_that("existing db connection passes assertion", {
  # testing this without setting a global variable is difficult;
  #. for now, we use something that will be present in the testthat execution environment
  # rather than trying to assign and test in the fixture itself
  thiscon <<- 1
  expect_silent(assert_db_connected(con="thiscon"))
})