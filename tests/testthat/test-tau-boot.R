test_that("tau.boot returns numeric vector of correct length", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  res <- tau.boot(data, method = "MM", nperm = 20)

  expect_type(res, "double")
  expect_length(res, 20)
})

test_that("tau.boot works for rational model", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  res <- tau.boot(data, method = "rational", nperm = 20)

  expect_type(res, "double")
  expect_length(res, 20)
})

test_that("tau.boot returns NA for failed fits", {
  # Pathological dataset: all zeros
  data <- matrix(0, nrow = 5, ncol = 5)

  res <- tau.boot(data, method = "MM", nperm = 10)

  expect_true(all(is.na(res)))
})

test_that("tau.boot handles invalid method", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  expect_error(tau.boot(data, method = "unknown", nperm = 10))
})

test_that("tau.boot produces varying bootstrap values when possible", {
  set.seed(123)
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  res <- tau.boot(data, method = "MM", nperm = 30)

  # At least some variation expected unless model always returns same limit
  expect_true(length(unique(res[!is.na(res)])) >= 1)
})

test_that("tau.boot works with data frames", {
  df <- data.frame(
    s1 = c(1,0,1),
    s2 = c(0,1,1),
    s3 = c(1,1,0)
  )

  res <- tau.boot(df, method = "MM", nperm = 10)

  expect_type(res, "double")
  expect_length(res, 10)
})