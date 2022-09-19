
test_that("Lengths of outputs should be the same", {
  increment <- seq(5, 101, by = 5)
  vector <- runif(106)
  output <- add_index_to_vector(increment, vector)
  expect_equal(length(output[[1]]), length(output[[2]]))
})
test_that("Throw error of increment is out-of-bounds of vec", {
  increment <- seq(5, 105, by = 5)
  vector <- runif(100)
  expect_error(add_index_to_vector(increment, vector))
})

test_that("Throw error if increment length less than 2", {
  increment <- seq(5, 105, by = 5)
  vector <- runif(100)
  expect_error(add_index_to_vector(increment, vector))
})
test_that("Throw error if not numeric increment", {
  increment <- as.character(seq(5, 100, by = 5))
  vector <- runif(100)
  expect_error(add_index_to_vector(increment, vector))
})
