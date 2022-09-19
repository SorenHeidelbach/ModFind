test_that("chunk_mapping doesn't return input data", {
  chunk_size <- 1e3
  dt <- dt_read_mapping()

  old_pos_chunk <- paste(dt$pos, dt$chunk, sep = "_")
  dt_overextending <- get_overextending_reads(dt, chunk_size)
  new_pos_chunk <- paste(
    dt_overextending$pos,
    dt_overextending$chunk,
    sep = "_")
  expect_true(!any(new_pos_chunk  %in% old_pos_chunk))
})

test_that("chunk_mapping should not create previously non-existing chunks", {
  chunk_size <- 1e3
  dt <- dt_read_mapping()
  old_chunks <- unique(dt$chunk)
  dt_overextending <- get_overextending_reads(dt, chunk_size)
  new_chunks <- unique(dt_overextending$chunk)
  expect_true(all(new_chunks  %in% old_chunks))
})

test_that("downsample actually removes rows", {
  chunk_size <- 1e3
  dt <- dt_read_mapping_unequal()
  dt_downsampled <- downsample(dt, chunk_size)
  
  expect_true(nrow(dt) > nrow(dt_downsampled))
})
