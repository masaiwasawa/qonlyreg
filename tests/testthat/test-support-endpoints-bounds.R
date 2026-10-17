test_that("qlm_bounds gives the same mean bounds from exact endpoints and identical support bounds", {

  qhat <- rbind(
    c(10, 20, 30),
    c(15, 25, 35)
  )
  tau <- c(0.25, 0.50, 0.75)

  from_endpoints <- qlm_bounds(
    qhat = qhat,
    tau = tau,
    y_lower = 0,
    y_upper = 100
  )

  from_bounds <- qlm_bounds(
    qhat = qhat,
    tau = tau,
    y_lower_bound = 0,
    y_upper_bound = 100
  )

  expect_equal(
    from_endpoints$m_lower,
    from_bounds$m_lower
  )
  expect_equal(
    from_endpoints$m_upper,
    from_bounds$m_upper
  )
})


test_that("exact endpoints dominate looser support bounds in mean bounding", {

  qhat <- matrix(c(10, 20, 30), nrow = 1)
  tau <- c(0.25, 0.50, 0.75)

  endpoint_only <- qlm_bounds(
    qhat = qhat,
    tau = tau,
    y_lower = 0,
    y_upper = 100
  )

  endpoint_and_loose_bounds <- qlm_bounds(
    qhat = qhat,
    tau = tau,
    y_lower = 0,
    y_upper = 100,
    y_lower_bound = -10,
    y_upper_bound = 120
  )

  expect_equal(
    endpoint_only$m_lower,
    endpoint_and_loose_bounds$m_lower
  )
  expect_equal(
    endpoint_only$m_upper,
    endpoint_and_loose_bounds$m_upper
  )
})


test_that("inconsistent endpoint and support-bound information is rejected", {

  qhat <- matrix(c(10, 20, 30), nrow = 1)
  tau <- c(0.25, 0.50, 0.75)

  expect_error(
    qlm_bounds(
      qhat = qhat,
      tau = tau,
      y_lower = 0,
      y_upper = 100,
      y_lower_bound = 5
    ),
    "y_lower"
  )

  expect_error(
    qlm_bounds(
      qhat = qhat,
      tau = tau,
      y_lower = 0,
      y_upper = 100,
      y_upper_bound = 90
    ),
    "y_upper"
  )
})


test_that("support bounds affect mean bounds but not Omega when no endpoint is supplied", {

  qhat <- rbind(
    c(10, 20, 30),
    c(20, 30, 40),
    c(30, 40, 50)
  )
  tau <- c(0.25, 0.50, 0.75)

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  n_g <- c(100, 120, 150)

  set.seed(1)
  region_100 <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower_bound = 0,
    y_upper_bound = 100,
    B = 100
  )

  set.seed(1)
  region_200 <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower_bound = 0,
    y_upper_bound = 200,
    B = 100
  )

  # Upper mean bounds must change.
  expect_false(
    isTRUE(all.equal(
      region_100$bounds$m_upper,
      region_200$bounds$m_upper
    ))
  )

  # Lower mean bounds do not depend on the upper support bound.
  expect_equal(
    region_100$bounds$m_lower,
    region_200$bounds$m_lower
  )

  # Crucial test:
  # support bounds must NOT enter finite-difference density estimation.
  expect_equal(
    region_100$Omega,
    region_200$Omega,
    tolerance = 1e-12
  )
})


test_that("a lower exact endpoint enters Omega whereas a lower support bound does not", {

  qhat <- rbind(
    c(10, 20, 30),
    c(20, 30, 40),
    c(30, 40, 50)
  )
  tau <- c(0.25, 0.50, 0.75)

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  n_g <- c(100, 120, 150)

  set.seed(2)
  endpoint_version <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower = 0,
    y_upper_bound = 100,
    B = 100
  )

  set.seed(2)
  bound_only_version <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower_bound = 0,
    y_upper_bound = 100,
    B = 100
  )

  # The same numerical lower value gives the same mean bounds.
  expect_equal(
    endpoint_version$bounds$m_lower,
    bound_only_version$bounds$m_lower
  )
  expect_equal(
    endpoint_version$bounds$m_upper,
    bound_only_version$bounds$m_upper
  )

  # But only the exact endpoint should enter finite differences.
  expect_false(
    isTRUE(all.equal(
      endpoint_version$Omega,
      bound_only_version$Omega
    ))
  )
})


test_that("changing an exact endpoint can change Omega", {

  qhat <- rbind(
    c(10, 20, 30),
    c(20, 30, 40),
    c(30, 40, 50)
  )
  tau <- c(0.25, 0.50, 0.75)

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  n_g <- c(100, 120, 150)

  set.seed(3)
  region_a <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower = 0,
    y_upper = 100,
    B = 100
  )

  set.seed(3)
  region_b <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower = 5,
    y_upper = 90,
    B = 100
  )

  expect_false(
    isTRUE(all.equal(
      region_a$Omega,
      region_b$Omega
    ))
  )
})


test_that("upper support bound changes bounds but not Omega when lower endpoint is exact", {

  qhat <- rbind(
    c(10, 20, 30),
    c(20, 30, 40),
    c(30, 40, 50)
  )
  tau <- c(0.25, 0.50, 0.75)

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  n_g <- c(100, 120, 150)

  set.seed(4)
  region_100 <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower = 0,
    y_upper_bound = 100,
    B = 100
  )

  set.seed(4)
  region_200 <- qlm_region(
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    y_lower = 0,
    y_upper_bound = 200,
    B = 100
  )

  # Identification changes.
  expect_false(
    isTRUE(all.equal(
      region_100$bounds$m_upper,
      region_200$bounds$m_upper
    ))
  )

  # Covariance plug-in should not change because U is only a bound.
  expect_equal(
    region_100$Omega,
    region_200$Omega,
    tolerance = 1e-12
  )
})


test_that("one disclosed quantile works with a lower exact endpoint", {

  qhat <- matrix(c(20, 30, 40), ncol = 1)
  tau <- 0.50

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  n_g <- c(100, 120, 150)

  expect_no_error(
    qlm_region(
      qhat = qhat,
      tau = tau,
      X = X,
      n_g = n_g,
      y_lower = 0,
      y_upper_bound = 100,
      B = 100
    )
  )
})


test_that("one disclosed quantile works with an upper exact endpoint", {

  qhat <- matrix(c(20, 30, 40), ncol = 1)
  tau <- 0.50

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  n_g <- c(100, 120, 150)

  expect_no_error(
    qlm_region(
      qhat = qhat,
      tau = tau,
      X = X,
      n_g = n_g,
      y_lower_bound = 0,
      y_upper = 100,
      B = 100
    )
  )
})


test_that("one disclosed quantile cannot use support bounds alone for finite-difference inference", {

  qhat <- matrix(c(20, 30, 40), ncol = 1)
  tau <- 0.50

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  n_g <- c(100, 120, 150)

  expect_error(
    qlm_region(
      qhat = qhat,
      tau = tau,
      X = X,
      n_g = n_g,
      y_lower_bound = 0,
      y_upper_bound = 100,
      B = 100
    ),
    "endpoint|finite-difference|quantile"
  )
})

test_that("legacy endpoint interface remains available", {

  qhat <- rbind(
    c(10, 20, 30),
    c(20, 30, 40)
  )
  tau <- c(0.25, 0.50, 0.75)

  expect_no_error(
    qlm_bounds(
      qhat = qhat,
      tau = tau,
      y_lower = 0,
      y_upper = 100
    )
  )

  X <- cbind(
    intercept = 1,
    x = c(0, 1)
  )

  expect_no_error(
    qlm_region(
      qhat = qhat,
      tau = tau,
      X = X,
      n_g = c(100, 120),
      y_lower = 0,
      y_upper = 100,
      B = 100
    )
  )
})


test_that("data.frame disclosure reads endpoint and support-bound columns separately", {

  disclosure <- data.frame(
    group = rep(1:2, each = 3),
    type = "quantile",
    tau = rep(c(0.25, 0.50, 0.75), 2),
    value = c(
      10, 20, 30,
      20, 30, 40
    ),
    n = rep(c(100, 120), each = 3),
    y_lower = c(
      rep(0, 3),
      rep(NA_real_, 3)
    ),
    y_lower_bound = c(
      rep(NA_real_, 3),
      rep(0, 3)
    ),
    y_upper_bound = 100
  )

  X <- cbind(
    intercept = 1,
    x = c(0, 1)
  )

  expect_no_error(
    qlm_region(
      disclosure = disclosure,
      X = X,
      B = 100
    )
  )
})


test_that("mixed disclosure works with support bounds for quantile groups", {

  disclosure <- data.frame(
    group = c(
      1, 1, 1,
      2,
      3, 3, 3
    ),
    type = c(
      "quantile", "quantile", "quantile",
      "mean",
      "quantile", "quantile", "quantile"
    ),
    tau = c(
      0.25, 0.50, 0.75,
      NA,
      0.25, 0.50, 0.75
    ),
    value = c(
      10, 20, 30,
      25,
      30, 40, 50
    ),
    n = c(
      100, 100, 100,
      120,
      150, 150, 150
    ),
    se = c(
      NA, NA, NA,
      1.5,
      NA, NA, NA
    )
  )

  X <- cbind(
    intercept = 1,
    x = c(0, 1, 2)
  )

  expect_no_error(
    qlm_region(
      disclosure = disclosure,
      X = X,
      y_lower_bound = 0,
      y_upper_bound = 100,
      B = 100
    )
  )
})


test_that("contrast bounds propagate support-bound arguments", {

  disclosure <- data.frame(
    study = rep(1, 6),
    arm = rep(c("case", "control"), each = 3),
    type = "quantile",
    tau = rep(c(0.25, 0.50, 0.75), 2),
    value = c(
      10, 20, 30,
      8, 18, 28
    )
  )

  expect_no_error(
    qlm_contrast_bounds(
      disclosure = disclosure,
      treatment = "case",
      control = "control",
      y_lower_bound = 0,
      y_upper_bound = 100
    )
  )
})


test_that("contrast region propagates endpoint and support-bound arguments", {

  disclosure <- data.frame(
    study = rep(1:2, each = 6),
    arm = rep(
      rep(c("case", "control"), each = 3),
      2
    ),
    type = "quantile",
    tau = rep(c(0.25, 0.50, 0.75), 4),
    value = c(
      10, 20, 30,
      8, 18, 28,
      15, 25, 35,
      12, 22, 32
    ),
    n = rep(c(100, 100, 100, 110, 110, 110), 2)
  )

  X <- cbind(
    intercept = 1,
    treatment = c(1, 0)
  )

  expect_no_error(
    qlm_contrast_region(
      disclosure = disclosure,
      X = X,
      treatment = "case",
      control = "control",
      y_lower = 0,
      y_upper_bound = 100,
      B = 100
    )
  )
})

test_that("support bounds are not partially matched as exact endpoints", {

  disclosure <- qonlyreg:::as_qlm_region_disclosure(
    qhat = matrix(c(10, 20, 30), nrow = 1),
    tau = c(0.25, 0.50, 0.75),
    n_g = 100,
    y_lower_bound = 0,
    y_upper_bound = 100
  )

  expect_null(disclosure[[1]][["y_lower"]])
  expect_null(disclosure[[1]][["y_upper"]])

  expect_equal(
    disclosure[[1]][["y_lower_bound"]],
    0
  )
  expect_equal(
    disclosure[[1]][["y_upper_bound"]],
    100
  )
})


test_that("qlm_bounds returns the correct quantile weights", {

  qhat <- matrix(
    c(10, 20, 30),
    nrow = 1
  )

  tau <- c(0.1, 0.4, 0.8)

  bounds <- qlm_bounds(
    qhat = qhat,
    tau = tau,
    y_lower_bound = 0,
    y_upper_bound = 40
  )

  delta_tau <- diff(c(0, tau, 1))

  expected_lower_weights <- delta_tau[-1L]
  expected_upper_weights <-
    delta_tau[-length(delta_tau)]

  expect_equal(
    bounds$weights_lower[[1]],
    expected_lower_weights
  )

  expect_equal(
    bounds$weights_upper[[1]],
    expected_upper_weights
  )

  # The returned weights must reconstruct the reported lower bound.
  expect_equal(
    tau[1] * 0 +
      sum(bounds$weights_lower[[1]] * qhat[1, ]),
    bounds$m_lower[1]
  )

  # The returned weights must reconstruct the reported upper bound.
  expect_equal(
    sum(bounds$weights_upper[[1]] * qhat[1, ]) +
      (1 - tau[length(tau)]) * 40,
    bounds$m_upper[1]
  )

  expect_equal(bounds$m_lower[1], 17)
  expect_equal(bounds$m_upper[1], 27)
})
