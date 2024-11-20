test_that("non-existant config location fails", {
  expect_error(check_config_location(config_file = "/home/turnip"))
})


test_that("weird config location raises warning", {
  expect_warning(check_config_location(config_file = "test_drugs_tasks.rds"))
})

test_that("non-existant warning", {
  config_file <- approved_locations[file.exists(approved_locations)][[1]]
  expect_silent(check_config_location(config_file = config_file))
})
