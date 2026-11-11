blocks <- list(
  A = simulate_lagged_series(n_times = 12, seed = 1)$X,
  B = simulate_lagged_series(n_times = 12, seed = 2)$X,
  C = simulate_lagged_series(n_times = 12, seed = 3)$X
)

test_that("all unordered pairs are analysed by default", {
  res <- lagged_protest_pairs(blocks, lags = 0:2, nperm = 49, seed = 1)
  expect_equal(nrow(res), 3)
  expect_equal(paste(res$x, res$y), c("A B", "A C", "B C"))
  expect_named(attr(res, "fits"), c("A -> B", "A -> C", "B -> C"))
  expect_equal(res$p.global.adjusted, p.adjust(res$p.global, "BH"))
  fits <- attr(res, "fits")
  expect_equal(res$Z, unname(vapply(fits, function(f) f$effect$Z, 1)))
  expect_equal(res$A, unname(vapply(fits, function(f) f$effect$A, 1)))
})

test_that("each row matches a direct call to lagged_protest()", {
  res <- lagged_protest_pairs(blocks, pairs = list(c("C", "A")), lags = 0:2,
                              nperm = 49, seed = 1)
  ref <- lagged_protest(blocks$C, blocks$A, lags = 0:2, bidirectional = TRUE,
                        nperm = 49, seed = 1)
  expect_equal(res$optimal_lag, ref$optimal_lag)
  expect_equal(res$m2, ref$m2)
  expect_equal(res$p.global, ref$p.global)
  expect_equal(res$n, 12 - abs(ref$optimal_lag))
  expect_equal(res$p.global.adjusted, res$p.global)
})

test_that("pairs can be given as a matrix and arguments are passed on", {
  res <- lagged_protest_pairs(blocks, pairs = cbind(c("A", "B"), c("B", "C")),
                              lags = 0:1, bidirectional = FALSE,
                              permutation = "cyclic")
  expect_equal(paste(res$x, res$y), c("A B", "B C"))
  expect_true(all(res$optimal_lag >= 0))
  expect_equal(attr(res, "fits")[[1]]$permutation, "cyclic")
})

test_that("the truncation message is emitted once", {
  long <- lapply(1:3, function(i)
    simulate_lagged_series(n_times = 14, seed = i)$X)
  names(long) <- c("A", "B", "C")
  expect_message(
    res <- lagged_protest_pairs(long, lags = 0:1, nperm = 39, period = 6,
                                seed = 1),
    "first 12 of 14 rows"
  )
  msgs <- testthat::capture_messages(
    lagged_protest_pairs(long, lags = 0:1, nperm = 39, period = 6, seed = 1)
  )
  expect_length(msgs, 1)
})

test_that("invalid inputs are rejected", {
  expect_error(lagged_protest_pairs(list(blocks$A, blocks$B)), "unique names")
  expect_error(lagged_protest_pairs(blocks["A"]), "at least two")
  expect_error(lagged_protest_pairs(blocks, pairs = list(c("A", "Z"))),
               "Unknown data sets")
  expect_error(lagged_protest_pairs(blocks, pairs = list(c("A", "A"))),
               "two different")
  expect_error(lagged_protest_pairs(blocks, pairs = matrix("A", 1, 3)),
               "two columns")
  expect_error(lagged_protest_pairs(blocks, pairs = list("A")),
               "length 2")
})
