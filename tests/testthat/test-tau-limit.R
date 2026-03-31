test_that("tau.limit returns correct structure for MM model", {
  data <- matrix(c(1,0,1,0, 1,1,0,0, 0,1,1,1), nrow = 3, byrow = TRUE)

  res <- tau.limit(data, method = "MM")

  expect_type(res, "list")
  expect_named(res, c("limit", "coef"))
  expect_type(res$limit, "double")
  expect_true(length(res$coef) == 2)
  expect_true(all(names(res$coef) %in% c("a", "b")))
})

test_that("tau.limit returns correct structure for rational model", {
  data <- matrix(c(1,0,1,0, 1,1,0,0, 0,1,1,1), nrow = 3, byrow = TRUE)

  res <- tau.limit(data, method = "rational")

  expect_type(res, "list")
  expect_named(res, c("limit", "coef"))
  expect_type(res$limit, "double")
  expect_true(length(res$coef) == 3)
  expect_true(all(names(res$coef) %in% c("a", "b", "c")))
})

test_that("tau.limit handles invalid method gracefully", {
  data <- matrix(c(1,0,1,0, 1,1,0,0, 0,1,1,1), nrow = 3, byrow = TRUE)

  expect_error(tau.limit(data, method = "unknown"))
})

test_that("tau.limit returns NA on fitting failure", {
  # Pathological dataset: all-zero rows
  data <- matrix(0, nrow = 5, ncol = 5)

  res1 <- tau.limit(data, method = "MM")
  res2 <- tau.limit(data, method = "rational")

  expect_true(is.na(res1$limit))
  expect_true(all(is.na(res1$coef)))

  expect_true(is.na(res2$limit))
  expect_true(all(is.na(res2$coef)))
})

test_that("tau.limit limit is >= observed richness when model converges", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  # MM model: should almost always converge on such data
  res <- tau.limit(data, method = "MM")
  expect_false(is.na(res$limit))
  expect_true(res$limit >= 3)

  # Rational model: may fail to converge; if it converges, limit >= observed richness
  res2 <- tau.limit(data, method = "rational")
  if (!is.na(res2$limit)) {
    expect_true(res2$limit >= 3)
  } else {
    succeed("Rational model did not converge on toy data, returning NA as designed.")
  }
})

test_that("tau.limit works with data frames as well as matrices", {
  df <- data.frame(
    s1 = c(1,0,1),
    s2 = c(0,1,1),
    s3 = c(1,1,0)
  )

  res <- tau.limit(df, method = "MM")
  expect_type(res$limit, "double")
})