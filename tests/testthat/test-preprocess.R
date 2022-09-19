# Testfile scope
# HDF5 file
hdf5_path <- "/home/ubuntu/test_data_modfind/PCR/signal_mappings.hdf5"
metainfo <- fread("/home/ubuntu/ModFind/data/metainfo_test_data.tsv")
h5_file <- hdf5r::H5File$new(hdf5_path, mode = "r")

# Slow loading and therefore only loaded once
dacs_vec <- h5_file[[glue::glue('/Batches/Batch_0/Dacs')]][]
ref_to_sig <-  h5_file[[glue::glue('/Batches/Batch_0/Ref_to_signal')]][]

test_that("add_signal returns list column", {
  metainfo_test <- copy(metainfo)
  add_signal(metainfo_test, dacs_vec = dacs_vec, ref_to_sig = ref_to_sig)
  expect_true(is.list(metainfo_test$current))
  expect_true(is.list(metainfo_test$current[[1]]))
  expect_true(is.list(metainfo_test$current[[2]]))
})

test_that("add_signal throws error if required columns are missing", {
  metainfo_test <- copy(metainfo)
  metainfo_test[, dacs_end := NULL]
  expect_error(
    add_signal(metainfo_test, dacs_vec = dacs_vec, ref_to_sig = ref_to_sig),
    "Following required columns are missing"
  )
})

test_that("add_signal throws error if there are duplicated read_ids", {
  metainfo_test <- rbind(metainfo, metainfo)
  expect_error(
    add_signal(metainfo_test, dacs_vec = dacs_vec, ref_to_sig = ref_to_sig),
    "Following rows have duplicate read_id"
  )
})

test_that("add_signal return same number of rows", {
  metainfo_test <- copy(metainfo)
  n_row_before <- nrow(metainfo_test)
  add_signal(metainfo_test, dacs_vec = dacs_vec, ref_to_sig = ref_to_sig)
  n_row_after <- nrow(metainfo_test)
  expect_equal(n_row_before, n_row_after)
})
