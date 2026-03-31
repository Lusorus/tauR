test_that("tau.plot runs without errors for MM model", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  expect_silent(
    tau.plot(data, method = "MM", new = TRUE)
  )
})

test_that("tau.plot runs without errors for rational model", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  expect_silent(
    tau.plot(data, method = "rational", new = TRUE)
  )
})

test_that("tau.plot works when tau.limit returns NA", {
  # Pathological dataset: all zeros → tau.limit() returns NA
  data <- matrix(0, nrow = 5, ncol = 5)

  expect_silent(
    tau.plot(data, method = "MM", new = TRUE)
  )
})

test_that("tau.plot adds to existing plot when new = FALSE", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  plot(1, type = "n")  # create empty plot

  expect_silent(
    tau.plot(data, method = "MM", new = FALSE)
  )
})

test_that("tau.plot accepts data frames", {
  df <- data.frame(
    s1 = c(1,0,1),
    s2 = c(0,1,1),
    s3 = c(1,1,0)
  )

  expect_silent(
    tau.plot(df, method = "MM", new = TRUE)
  )
})

test_that("tau.plot respects xlim when provided", {
  data <- matrix(c(1,0,1,0,
                   1,1,0,0,
                   0,1,1,1),
                 nrow = 3, byrow = TRUE)

  expect_silent(
    tau.plot(data, method = "MM", new = TRUE, xlim = c(1, 10))
  )
})