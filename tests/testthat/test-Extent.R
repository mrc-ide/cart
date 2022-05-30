test_that("Get extent works", {
  r <- terra::rast()
  extent_in <- c(0, 2.5, 0, 1.5)
  terra::ext(r) <- extent_in
  expect_identical(get_extent(r), matrix(extent_in, ncol = 2, byrow = TRUE))
})
