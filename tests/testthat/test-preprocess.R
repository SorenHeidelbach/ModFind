# Testfile scope
# HDF5 file
hdf5_path <- "/home/ubuntu/test_data_modfind/PCR/signal_mappings.hdf5"
metainfo <- fread("/home/ubuntu/nanomotif/data/metainfo_test_data.tsv")
h5_file <- hdf5r::H5File$new(hdf5_path, mode = "r")

test_that("add_signal returns list column", {
  metainfo_test <- copy(metainfo)
  add_signal(metainfo_test, h5_file, "Batch_0")
  expect_true(is.list(metainfo_test$current))
  expect_true(is.list(metainfo_test$current[[1]]))
  expect_true(is.list(metainfo_test$current[[2]]))
})

test_that("add_signal throws error if required columns are missing", {
  metainfo_test <- copy(metainfo)
  metainfo_test[, dacs_end := NULL]
  expect_error(
    add_signal(metainfo_test, h5_file, "Batch_0"),
    "Following required columns are missing"
  )
})

test_that("add_signal throws error if there are duplicated read_ids", {
  metainfo_test <- rbind(metainfo, metainfo)
  expect_error(
    add_signal(metainfo_test, h5_file, "Batch_0"),
    "Following rows have duplicate read_id"
  )
})

test_that("add_signal return same number of rows", {
  metainfo_test <- copy(metainfo)
  n_row_before <- nrow(metainfo_test)
  add_signal(metainfo_test, h5_file, "Batch_0")
  n_row_after <- nrow(metainfo_test)
  expect_equal(n_row_before, n_row_after)
})

