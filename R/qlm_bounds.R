#' Compute mean bounds from disclosed quantiles
#'
#' Compute sharp lower and upper bounds for group conditional means
#' using disclosed group-level quantiles and known support endpoints.
#'
#' @param qhat A G by K matrix of disclosed sample quantiles. Each row
#'   corresponds to a group and each column corresponds to a quantile
#'   level in `tau`. Quantiles must be ordered increasingly within each row.
#' @param tau A numeric vector of length K containing quantile levels
#'   in (0, 1), sorted in increasing order.
#' @param y_lower A numeric vector of length G, or a scalar, giving known
#'   lower support endpoints for each group.
#' @param y_upper A numeric vector of length G, or a scalar, giving known
#'   upper support endpoints for each group.
#'
#' @return
#' An object of class `"qlm_bounds"` containing:
#' \describe{
#'   \item{m_lower}{Estimated lower bounds for group conditional means.}
#'   \item{m_upper}{Estimated upper bounds for group conditional means.}
#'   \item{qhat}{Input matrix of disclosed sample quantiles.}
#'   \item{tau}{Input quantile levels.}
#'   \item{y_lower}{Lower support endpoints.}
#'   \item{y_upper}{Upper support endpoints.}
#'   \item{weights_lower}{Weights on disclosed quantiles for lower bounds.}
#'   \item{weights_upper}{Weights on disclosed quantiles for upper bounds.}
#' }
#'
#' @export
qlm_bounds <- function(qhat, tau, y_lower, y_upper) {

  qhat <- as.matrix(qhat)
  tau <- as.numeric(tau)

  G <- nrow(qhat)
  K <- ncol(qhat)

  if (length(tau) != K) {
    stop("length(tau) must equal ncol(qhat).")
  }

  if (any(tau <= 0 | tau >= 1)) {
    stop("All elements of tau must lie in (0, 1).")
  }

  if (is.unsorted(tau)) {
    stop("tau must be sorted in increasing order.")
  }

  if (any(apply(qhat, 1, is.unsorted))) {
    stop("Each row of qhat must be sorted in increasing order.")
  }

  if (length(y_lower) == 1) {
    y_lower <- rep(y_lower, G)
  }

  if (length(y_upper) == 1) {
    y_upper <- rep(y_upper, G)
  }

  y_lower <- as.numeric(y_lower)
  y_upper <- as.numeric(y_upper)

  if (length(y_lower) != G) {
    stop("y_lower must be a scalar or a vector of length nrow(qhat).")
  }

  if (length(y_upper) != G) {
    stop("y_upper must be a scalar or a vector of length nrow(qhat).")
  }

  if (any(y_lower > y_upper)) {
    stop("Each element of y_lower must be less than or equal to y_upper.")
  }

  if (any(qhat[, 1] < y_lower)) {
    stop("The first disclosed quantile must be greater than or equal to y_lower.")
  }

  if (any(qhat[, K] > y_upper)) {
    stop("The last disclosed quantile must be less than or equal to y_upper.")
  }

  tau_full <- c(0, tau, 1)
  delta_tau <- diff(tau_full)

  weights_lower <- delta_tau[-1]
  weights_upper <- delta_tau[-length(delta_tau)]

  m_lower <- delta_tau[1] * y_lower +
    as.numeric(qhat %*% weights_lower)

  m_upper <- as.numeric(qhat %*% weights_upper) +
    delta_tau[length(delta_tau)] * y_upper

  out <- list(
    m_lower = m_lower,
    m_upper = m_upper,
    qhat = qhat,
    tau = tau,
    y_lower = y_lower,
    y_upper = y_upper,
    weights_lower = weights_lower,
    weights_upper = weights_upper
  )

  class(out) <- "qlm_bounds"
  out
}


#' @export
print.qlm_bounds <- function(x, ...) {
  cat("Quantile-only mean bounds\n\n")
  out <- data.frame(
    m_lower = x$m_lower,
    m_upper = x$m_upper
  )
  print(out)
  invisible(x)
}
