qlm_omega_bounds <- function(qhat, tau, n_g, density_floor = 1e-8) {

  qhat <- as.matrix(qhat)
  tau <- as.numeric(tau)
  n_g <- as.numeric(n_g)

  G <- nrow(qhat)

  tau_full <- c(0, tau, 1)
  delta_tau <- diff(tau_full)

  aL <- delta_tau[-1]
  aU <- delta_tau[-length(delta_tau)]

  Omega <- matrix(0, nrow = 2 * G, ncol = 2 * G)

  for (g in seq_len(G)) {

    fhat <- qlm_density_at_quantiles(
      qhat = qhat[g, ],
      tau = tau,
      density_floor = density_floor
    )

    cov_tau <- outer(tau, tau, function(a, b) pmin(a, b) - a * b)
    denom <- outer(fhat, fhat, "*")
    Sigma_g <- cov_tau / denom

    Omega_g <- (sum(n_g) / n_g[g]) *
      matrix(
        c(
          as.numeric(t(aL) %*% Sigma_g %*% aL),
          as.numeric(t(aL) %*% Sigma_g %*% aU),
          as.numeric(t(aU) %*% Sigma_g %*% aL),
          as.numeric(t(aU) %*% Sigma_g %*% aU)
        ),
        nrow = 2,
        byrow = TRUE
      )

    idx <- c(2 * g - 1, 2 * g)
    Omega[idx, idx] <- Omega_g
  }

  Omega
}


#' Simultaneous confidence region from quantile-only partial identification
#'
#' Construct a simultaneous confidence region for regression coefficients
#' using disclosed group-level quantiles.
#'
#' @param qhat A G by K matrix of disclosed sample quantiles.
#'   Each row corresponds to a group and each column corresponds
#'   to a quantile level in `tau`. Quantiles must be ordered
#'   increasingly within each row.
#' @param tau A numeric vector of length K containing quantile
#'   levels in (0, 1), sorted in increasing order.
#' @param X A G by p matrix of group-level covariates.
#' @param n_g A numeric vector of length G containing group sizes.
#' @param y_lower A numeric vector of length G, or a scalar, giving known
#'   lower support endpoints for each group.
#' @param y_upper A numeric vector of length G, or a scalar, giving known
#'   upper support endpoints for each group.
#' @param level Confidence level. Default is 0.95.
#' @param B Number of Gaussian simulation draws used to compute the
#'   simultaneous critical value.
#' @param density_floor A positive lower bound used to avoid
#'   numerical instability in density estimation.
#'
#' @return
#' An object of class `"qlm_region"` containing:
#' \describe{
#'   \item{bounds}{Estimated group mean bounds.}
#'   \item{Omega}{Estimated covariance matrix of the estimated mean bounds.}
#'   \item{critical_value}{Simulated simultaneous critical value.}
#'   \item{level}{Confidence level.}
#'   \item{lower_adjusted}{Adjusted lower constraints.}
#'   \item{upper_adjusted}{Adjusted upper constraints.}
#'   \item{const_mat}{Constraint matrix defining the simultaneous confidence region.}
#'   \item{const_dir}{Constraint directions.}
#'   \item{const_rhs}{Right-hand side vector of the constraints.}
#'   \item{call}{Matched function call.}
#' }
#'
#' @examples
#' set.seed(123)
#'
#' N <- 200
#' G <- 5
#' tau <- c(0.25, 0.5, 0.75)
#'
#' X_grid <- seq(-1, 1, length.out = G)
#' g <- rep(1:G, length.out = N)
#' Xg <- X_grid[g]
#' Y <- 1 + 2 * Xg + rnorm(N)
#'
#' qhat_mat <- t(sapply(seq_len(G), function(j) {
#'   stats::quantile(Y[g == j], probs = tau, type = 1)
#' }))
#'
#' X <- stats::model.matrix(~ X_grid)
#' n_g <- as.numeric(table(g))
#'
#' region <- qlm_region(
#'   qhat = qhat_mat,
#'   tau = tau,
#'   X = X,
#'   n_g = n_g,
#'   y_lower = -10,
#'   y_upper = 10,
#'   B = 100
#' )
#'
#' confint(region)
#'
#' @export
qlm_region <- function(qhat, tau, X, n_g,
                       y_lower, y_upper,
                       level = 0.95,
                       B = 10000,
                       density_floor = 1e-8) {

  qhat <- as.matrix(qhat)
  X <- as.matrix(X)
  tau <- as.numeric(tau)
  n_g <- as.numeric(n_g)

  G <- nrow(qhat)
  N <- sum(n_g)

  if (nrow(X) != G) {
    stop("nrow(X) must equal nrow(qhat).")
  }

  bounds <- qlm_bounds(
    qhat = qhat,
    tau = tau,
    y_lower = y_lower,
    y_upper = y_upper
  )

  Omega <- qlm_omega_bounds(
    qhat = qhat,
    tau = tau,
    n_g = n_g,
    density_floor = density_floor
  )

  omega <- sqrt(diag(Omega))

  if (any(omega <= 0)) {
    stop(
      "Estimated asymptotic variances must be strictly positive. ",
      "This may occur when disclosed quantiles are nearly identical."
    )
  }

  Z <- MASS::mvrnorm(
    n = B,
    mu = rep(0, 2 * G),
    Sigma = Omega
  )

  Tsim <- apply(abs(sweep(Z, 2, omega, "/")), 1, max)
  qhat_crit <- as.numeric(stats::quantile(Tsim, probs = level, type = 8))

  adj <- qhat_crit * omega / sqrt(N)

  lower_adj <- bounds$m_lower - adj[seq(1, 2 * G, by = 2)]
  upper_adj <- bounds$m_upper + adj[seq(2, 2 * G, by = 2)]

  const_mat <- rbind(X, X)
  const_dir <- c(rep(">=", G), rep("<=", G))
  const_rhs <- c(lower_adj, upper_adj)

  out <- list(
    bounds = bounds,
    Omega = Omega,
    critical_value = qhat_crit,
    level = level,
    lower_adjusted = lower_adj,
    upper_adjusted = upper_adj,
    const_mat = const_mat,
    const_dir = const_dir,
    const_rhs = const_rhs,
    X = X,
    n_g = n_g,
    call = match.call()
  )

  class(out) <- "qlm_region"
  out
}


#' @export
print.qlm_region <- function(x, ...) {
  cat("Quantile-only simultaneous confidence region\n\n")
  cat("Confidence level:", x$level, "\n")
  cat("Critical value:", x$critical_value, "\n\n")
  cat("Linear constraints:\n")
  cat("  const_mat %*% beta  const_dir  const_rhs\n\n")
  cat("Number of constraints:", nrow(x$const_mat), "\n")
  cat("Number of coefficients:", ncol(x$const_mat), "\n")
  invisible(x)
}


#' @export
summary.qlm_region <- function(object, ...) {

  out <- list(
    call = object$call,
    level = object$level,
    critical_value = object$critical_value,
    lower_adjusted = object$lower_adjusted,
    upper_adjusted = object$upper_adjusted,
    const_mat = object$const_mat,
    const_dir = object$const_dir,
    const_rhs = object$const_rhs
  )

  class(out) <- "summary.qlm_region"
  out
}


#' @export
print.summary.qlm_region <- function(x, ...) {
  cat("Quantile-only simultaneous confidence region\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nConfidence level:", x$level, "\n")
  cat("Critical value:", x$critical_value, "\n\n")

  cat("Adjusted group mean constraints:\n")
  print(
    data.frame(
      lower = x$lower_adjusted,
      upper = x$upper_adjusted
    )
  )

  invisible(x)
}


#' Confidence intervals from a qlm region
#'
#' Compute coefficient-wise confidence intervals by projecting a
#' simultaneous confidence region onto each coordinate axis.
#'
#' @param object An object of class `"qlm_region"`.
#' @param parm Currently ignored.
#' @param level Currently ignored. The confidence level is determined
#'   when constructing the `"qlm_region"` object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return
#' An object of class `"qlm_confint"` containing:
#' \describe{
#'   \item{confint}{Coefficient-wise confidence intervals.}
#'   \item{level}{Confidence level.}
#'   \item{critical_value}{Simulated simultaneous critical value.}
#'   \item{call}{Matched function call of the original region object.}
#' }
#'
#' @export
confint.qlm_region <- function(object, parm = NULL, level = NULL, ...) {

  const_mat <- object$const_mat
  const_dir <- object$const_dir
  const_rhs <- object$const_rhs

  p <- ncol(const_mat)

  const_mat_free <- cbind(const_mat, -const_mat)

  ci <- matrix(NA_real_, nrow = p, ncol = 2)
  colnames(ci) <- c("lower", "upper")
  rownames(ci) <- colnames(const_mat)

  if (is.null(rownames(ci))) {
    rownames(ci) <- paste0("beta", seq_len(p))
  }

  for (j in seq_len(p)) {

    obj <- rep(0, 2 * p)
    obj[j] <- 1
    obj[p + j] <- -1

    sol_min <- lpSolve::lp(
      direction = "min",
      objective.in = obj,
      const.mat = const_mat_free,
      const.dir = const_dir,
      const.rhs = const_rhs
    )

    sol_max <- lpSolve::lp(
      direction = "max",
      objective.in = obj,
      const.mat = const_mat_free,
      const.dir = const_dir,
      const.rhs = const_rhs
    )

    if (sol_min$status != 0 || sol_max$status != 0) {
      warning("Linear program did not converge for coefficient ", j)
    } else {
      ci[j, 1] <- sol_min$objval
      ci[j, 2] <- sol_max$objval
    }
  }

  out <- list(
    confint = ci,
    level = object$level,
    critical_value = object$critical_value,
    call = object$call
  )

  class(out) <- "qlm_confint"
  out
}


#' @export
print.qlm_confint <- function(x, ...) {
  cat("Quantile-only coefficient-wise confidence intervals\n\n")
  cat("Confidence level:", x$level, "\n")
  cat("Critical value:", x$critical_value, "\n\n")
  print(x$confint)
  invisible(x)
}


#' @export
summary.qlm_confint <- function(object, ...) {

  out <- list(
    call = object$call,
    level = object$level,
    critical_value = object$critical_value,
    confint = object$confint
  )

  class(out) <- "summary.qlm_confint"
  out
}


#' @export
print.summary.qlm_confint <- function(x, ...) {
  cat("Quantile-only coefficient-wise confidence intervals\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nConfidence level:", x$level, "\n")
  cat("Critical value:", x$critical_value, "\n\n")

  cat("Coefficient-wise confidence intervals:\n")
  print(x$confint)

  invisible(x)
}


#' Plot coefficient-wise confidence intervals
#'
#' Plot coefficient-wise confidence intervals obtained by projecting
#' a quantile-only simultaneous confidence region.
#'
#' @param x An object of class `"qlm_confint"`.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A base R plot is produced. The function returns `x` invisibly.
#'
#' @export
plot.qlm_confint <- function(x, ...) {

  ci <- x$confint
  coef_names <- rownames(ci)

  if (is.null(coef_names)) {
    coef_names <- paste0("beta", seq_len(nrow(ci)))
  }

  y_pos <- seq_len(nrow(ci))
  xlim <- range(ci, finite = TRUE)

  graphics::plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, nrow(ci) + 0.5),
    yaxt = "n",
    xlab = "Coefficient value",
    ylab = "",
    main = paste0("Coefficient-wise confidence intervals, level = ", x$level)
  )

  graphics::axis(
    side = 2,
    at = y_pos,
    labels = coef_names,
    las = 1
  )

  graphics::segments(
    x0 = ci[, 1],
    y0 = y_pos,
    x1 = ci[, 2],
    y1 = y_pos,
    lwd = 2
  )

  graphics::points(
    x = rowMeans(ci),
    y = y_pos,
    pch = 19
  )

  graphics::abline(v = 0, lty = 2)

  invisible(x)
}
