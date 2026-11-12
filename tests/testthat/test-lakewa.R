test_that("lakewa_blocks() returns three aligned blocks without gaps", {
  skip_if_not_installed("MARSS")
  lw <- lakewa_blocks()
  expect_named(lw, c("Environment", "Phytoplankton", "Zooplankton"))
  expect_equal(vapply(lw, ncol, integer(1)),
               c(Environment = 3L, Phytoplankton = 5L, Zooplankton = 7L))
  expect_true(all(vapply(lw, nrow, integer(1)) == 396L))
  expect_identical(rownames(lw$Environment)[c(1, 396)],
                   c("1962-01", "1994-12"))
  expect_false(anyNA(unlist(lw)))
})

test_that("lakewa_blocks() matches the raw MARSS data", {
  skip_if_not_installed("MARSS")
  e <- new.env()
  utils::data("lakeWAplankton", package = "MARSS", envir = e)
  raw <- as.data.frame(e$lakeWAplanktonRaw)
  lw <- lakewa_blocks(interpolate = FALSE)
  expect_equal(unname(lw$Environment[, "Temp"]), raw$Temp)
  expect_equal(unname(lw$Zooplankton[, "Daphnia"]), log1p(raw$Daphnia))
  expect_true(anyNA(lw$Zooplankton))
})

test_that("lakewa_blocks() subsets years and validates input", {
  skip_if_not_installed("MARSS")
  lw <- lakewa_blocks(from = 1970, to = 1979)
  expect_equal(nrow(lw$Phytoplankton), 120L)
  expect_identical(rownames(lw$Phytoplankton)[1], "1970-01")
  expect_error(lakewa_blocks(from = 1950))
  expect_error(lakewa_blocks(from = 1990, to = 1980))
  expect_error(lakewa_blocks(interpolate = NA))
})

test_that("lakewa_blocks() explains how to get MARSS when it is missing", {
  local_mocked_bindings(has_marss = function() FALSE)
  expect_error(lakewa_blocks(), "needs the MARSS package")
})

test_that("fill_linear() interpolates inside and holds the ends", {
  expect_equal(fill_linear(c(NA, 1, NA, 3, NA)), c(1, 1, 2, 3, 3))
  expect_error(fill_linear(c(NA, 2, NA)), "Too few observed values")
})
