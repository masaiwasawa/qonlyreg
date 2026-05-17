#' Quantile-only linear model
#'
#' Estimate a linear regression model using group-level disclosed quantiles.
#'
#' @param qhat A G by K matrix of disclosed sample quantiles.
#'   Each row corresponds to a group and each column corresponds
#'   to a quantile level in `tau`. Quantiles must be ordered
#'   increasingly within each row.
#' @param tau A numeric vector of length K containing quantile
#'   levels in (0, 1), sorted in increasing order.
#' @param X A G by p matrix of group-level covariates.
#' @param n_g A numeric vector of length G containing
#'   group sizes.
#' @param density_floor A positive lower bound used to avoid
#'   numerical instability in density estimation.
#' @return
#' An object of class `"qlm"` containing:
#' \describe{
#'   \item{coefficients}{Estimated regression coefficients.}
#'   \item{vcov}{Estimated covariance matrix of the estimator.}
#'   \item{se}{Estimated standard errors.}
#'   \item{m_hat}{Estimated group conditional means.}
#'   \item{var_m_hat}{Estimated variances of the group mean estimators.}
#'   \item{weights}{Integration weights used for quantile aggregation.}
#'   \item{fhat}{Estimated densities at disclosed quantiles.}
#' }
#'
#' @export
qlm <- function(qhat, tau, X, n_g, density_floor = 1e-8) {

  qhat <- as.matrix(qhat)
  X <- as.matrix(X)
  tau <- as.numeric(tau)
  n_g <- as.numeric(n_g)

  G <- nrow(qhat)
  K <- ncol(qhat)

  if (length(tau) != K) {
    stop("length(tau) must equal ncol(qhat).")
  }
  if (nrow(X) != G) {
    stop("nrow(X) must equal nrow(qhat).")
  }
  if (length(n_g) != G) {
    stop("length(n_g) must equal nrow(qhat).")
  }
  if (any(tau <= 0 | tau >= 1)) {
    stop("All elements of tau must lie in (0, 1).")
  }
  if (is.unsorted(tau)) {
    stop("tau must be sorted in increasing order.")
  }
  if (any(n_g <= 0)) {
    stop("All group sizes n_g must be positive.")
  }
  if (any(apply(qhat, 1, is.unsorted))) {
    stop("Each row of qhat must be sorted in increasing order.")
  }

  alpha <- qlm_weights(tau)

  m_hat <- as.numeric(qhat %*% alpha)

  var_m_hat <- numeric(G)
  fhat <- matrix(NA_real_, nrow = G, ncol = K)

  for (g in seq_len(G)) {
    fhat[g, ] <- qlm_density_at_quantiles(
      qhat = qhat[g, ],
      tau = tau,
      density_floor = density_floor
    )

    var_m_hat[g] <- qlm_var_mean_bound(
      fhat = fhat[g, ],
      tau = tau,
      alpha = alpha,
      n_g = n_g[g]
    )
  }

  W <- diag(n_g, nrow = G)
  A <- crossprod(X, W %*% X)

  beta_hat <- solve(A, crossprod(X, W %*% m_hat))

  meat <- matrix(0, nrow = ncol(X), ncol = ncol(X))

  for (g in seq_len(G)) {
    xg <- matrix(X[g, ], ncol = 1)
    meat <- meat + n_g[g]^2 * var_m_hat[g] * (xg %*% t(xg))
  }

  vcov_beta <- solve(A, meat %*% solve(A))
  se_beta <- sqrt(diag(vcov_beta))

  out <- list(
    call = match.call(),
    coefficients = as.numeric(beta_hat),
    vcov = vcov_beta,
    se = se_beta,
    m_hat = m_hat,
    var_m_hat = var_m_hat,
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = n_g,
    weights = alpha,
    fhat = fhat
  )

  class(out) <- "qlm"
  out
}


#' @export
coef.qlm <- function(object, ...) {
  object$coefficients
}


#' @export
vcov.qlm <- function(object, ...) {
  object$vcov
}


#' @export
print.qlm <- function(x, ...) {
  cat("Quantile-only linear model\n\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  invisible(x)
}


#' @export
summary.qlm <- function(object, ...) {
  est <- object$coefficients
  se <- object$se
  z <- est / se
  pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

  coef_table <- cbind(
    Estimate = est,
    Std.Error = se,
    z.value = z,
    p.value = pval
  )

  out <- list(
    call = object$call,
    coefficients = coef_table
  )

  class(out) <- "summary.qlm"
  out
}


#' @export
print.summary.qlm <- function(x, ...) {
  cat("Quantile-only linear model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}


qlm_weights <- function(tau) {
  K <- length(tau)
  alpha <- numeric(K)

  if (K == 1) {
    alpha[1] <- 1
    return(alpha)
  }

  alpha[1] <- tau[1] + (tau[2] - tau[1]) / 2
  alpha[K] <- (tau[K] - tau[K - 1]) / 2 + (1 - tau[K])

  if (K >= 3) {
    for (k in 2:(K - 1)) {
      alpha[k] <- (tau[k] - tau[k - 1]) / 2 + (tau[k + 1] - tau[k]) / 2
    }
  }
  alpha
}


qlm_density_at_quantiles <- function(qhat, tau, density_floor = 1e-8) {
  K <- length(tau)

  if (K < 2) {
    stop("At least two disclosed quantiles are needed to estimate densities.")
  }

  qprime_hat <- numeric(K)

  qprime_hat[1] <- (qhat[2] - qhat[1]) / (tau[2] - tau[1])
  qprime_hat[K] <- (qhat[K] - qhat[K - 1]) / (tau[K] - tau[K - 1])

  if (K >= 3) {
    for (k in 2:(K - 1)) {
      qprime_hat[k] <- (qhat[k + 1] - qhat[k - 1]) /
        (tau[k + 1] - tau[k - 1])
    }
  }

  qprime_hat <- pmax(qprime_hat, density_floor)

  1 / qprime_hat
}


qlm_var_mean_bound <- function(fhat, tau, alpha, n_g) {
  cov_tau <- outer(tau, tau, function(a, b) pmin(a, b) - a * b)
  denom <- outer(fhat, fhat, "*")
  var_q <- cov_tau / denom / n_g

  as.numeric(t(alpha) %*% var_q %*% alpha)
}
