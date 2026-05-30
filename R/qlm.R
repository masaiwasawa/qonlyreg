#' Quantile/mean disclosure linear model
#'
#' Estimate a linear regression model using group-level disclosed quantiles,
#' observed means, or mixed disclosure.
#'
#' The function first constructs an estimated group conditional mean for each
#' group. For quantile disclosures, the group mean is approximated by numerical
#' integration of the disclosed quantile function. For mean disclosures, the
#' reported mean is used directly. The regression coefficient is then estimated
#' by weighted least squares using group sizes as weights.
#'
#' @param qhat A G by K matrix of disclosed sample quantiles. Each row
#'   corresponds to a group and each column corresponds to a quantile level
#'   in `tau`. This is the legacy common-grid quantile-only interface.
#'   Ignored if `disclosure` is supplied.
#' @param tau A numeric vector of length K containing quantile levels in
#'   (0, 1), sorted in increasing order. This is the legacy common-grid
#'   quantile-only interface. Ignored if `disclosure` is supplied.
#' @param X A G by p matrix of group-level covariates. The rows of `X` must
#'   correspond to the groups represented by `qhat` or `disclosure`, in the
#'   same order.
#' @param n_g A numeric vector of length G containing group sizes. This is
#'   the legacy common-grid quantile-only interface. Ignored if `disclosure`
#'   is supplied.
#' @param disclosure Optional group-level disclosure object. This may be
#'   either a list or a data.frame. If supplied, this takes precedence over
#'   `qhat`, `tau`, and `n_g`.
#'
#'   A data.frame input should be in long format. It must contain a group
#'   identifier column, a `type` column, and an `n` column. The group
#'   identifier column is detected from one of `group`, `group_id`, `study`,
#'   `study_id`, or `id`. The `type` column indicates whether each row
#'   reports a quantile or a mean. Supported values are `quantile`, `q`,
#'   `quantiles`, `mean`, and `m`.
#'
#'   Quantile rows require a `tau` column and one disclosed-value column.
#'   The disclosed-value column may be named `value`, `qhat`, or `estimate`.
#'
#'   Mean rows require either a `mean` column or one disclosed-value column
#'   named `value`, `qhat`, or `estimate`.
#'   For inference, at least one of `var`, `se`, or `sd` may additionally
#'   be supplied. If none of these quantities is provided for a mean
#'   disclosure, the coefficient covariance matrix and standard errors
#'   cannot be computed.
#'    Here `var` refers to the
#'   variance of the reported mean estimator, not the variance
#'   of the underlying observations.
#'
#'   A list input should contain one element per group. Each element should
#'   be a list with `type = "quantile"`, `type = "mean"`, or
#'   `type = "mean_quantile"`. Quantile disclosures require `tau`, `qhat`,
#'   and `n`. Mean disclosures require `mean` and `n`. For mean
#'   disclosures, at least one of `var`, `se`, or `sd` may additionally
#'   be supplied.
#' @param density_floor A positive lower bound used to avoid numerical
#'   instability in density estimation.
#'
#' @return
#' An object of class `"qlm"` containing:
#' \describe{
#'   \item{call}{The matched function call.}
#'   \item{coefficients}{Estimated regression coefficients.}
#'   \item{vcov}{Estimated covariance matrix of the estimator. If a mean
#'     disclosure is supplied without `var`, `se`, or `sd`, this matrix
#'     is filled with `NA`.}
#'   \item{se}{Estimated standard errors. If the covariance matrix cannot be
#'     computed, these are `NA`.}
#'   \item{m_hat}{Estimated group conditional means.}
#'   \item{var_m_hat}{Estimated variances of the group mean estimators.}
#'   \item{disclosure}{Standardized group-level disclosure objects.}
#'   \item{X}{Group-level covariate matrix.}
#'   \item{n_g}{Group sizes.}
#'   \item{weights}{A list of integration weights used for quantile aggregation.}
#'   \item{fhat}{A list of estimated densities at disclosed quantiles.}
#' }
#'
#' @examples
#' ## Legacy common-grid quantile-only interface
#' qhat <- rbind(
#'   c(10, 15, 20),
#'   c(20, 25, 30),
#'   c(30, 35, 40)
#' )
#'
#' tau <- c(0.25, 0.50, 0.75)
#'
#' X <- cbind(
#'   intercept = 1,
#'   x = c(0, 1, 2)
#' )
#'
#' n_g <- c(100, 120, 150)
#'
#' fit <- qlm(
#'   qhat = qhat,
#'   tau = tau,
#'   X = X,
#'   n_g = n_g
#' )
#'
#' fit
#' summary(fit)
#'
#' ## Disclosure data.frame interface with mixed disclosure
#' disclosure <- data.frame(
#'   group = c(
#'     1, 1, 1,
#'     2,
#'     3, 3, 3
#'   ),
#'   type = c(
#'     "quantile", "quantile", "quantile",
#'     "mean",
#'     "quantile", "quantile", "quantile"
#'   ),
#'   tau = c(
#'     0.25, 0.50, 0.75,
#'     NA,
#'     0.25, 0.50, 0.75
#'   ),
#'   value = c(
#'     10, 15, 20,
#'     25,
#'     30, 35, 40
#'   ),
#'   n = c(
#'     100, 100, 100,
#'     120,
#'     150, 150, 150
#'   ),
#'   se = c(
#'     NA, NA, NA,
#'     1.5,
#'     NA, NA, NA
#'   )
#' )
#'
#' fit_mixed <- qlm(
#'   disclosure = disclosure,
#'   X = X
#' )
#'
#' fit_mixed
#' summary(fit_mixed)
#'
#' @export
qlm <- function(qhat = NULL,
                tau = NULL,
                X,
                n_g = NULL,
                disclosure = NULL,
                density_floor = 1e-8) {

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }
  if (any(!is.finite(X))) {
    stop("X contains non-finite values.")
  }
  if (!is.numeric(density_floor) ||
      length(density_floor) != 1 ||
      !is.finite(density_floor) ||
      density_floor <= 0) {
    stop("density_floor must be a positive finite scalar.")
  }

  disclosure <- as_qlm_disclosure(
    disclosure = disclosure,
    qhat = qhat,
    tau = tau,
    n_g = n_g
  )

  G <- length(disclosure)
  p <- ncol(X)

  if (nrow(X) != G) {
    stop("nrow(X) must equal the number of groups in disclosure.")
  }
  if (G <= p) {
    warning(
      "The number of groups is not larger than the number of regressors. ",
      "The regression may be exactly determined or unstable."
    )
  }

  n_vec <- vapply(disclosure, function(d) d$n, numeric(1))

  qlm_warn_disclosure_structure(disclosure)

  m_hat <- numeric(G)
  var_m_hat <- rep(NA_real_, G)
  weights <- vector("list", G)
  fhat <- vector("list", G)

  for (g in seq_len(G)) {
    d <- disclosure[[g]]

    if (d$type %in% c("quantile", "mean_quantile")) {
      if (length(d$tau) < 2) {
        stop(
          "qlm() requires at least two disclosed quantiles for quantile-based ",
          "point inference. For partial identification with a single disclosed ",
          "quantile, use qlm_region() with known support bounds."
        )
      }

      alpha_g <- qlm_weights(d$tau)

      fhat_g <- qlm_density_at_quantiles(
        qhat = d$qhat,
        tau = d$tau,
        density_floor = density_floor
      )

      qlm_warn_density(
        fhat = fhat_g,
        qhat = d$qhat,
        tau = d$tau,
        group = g,
        density_floor = density_floor
      )

      m_quantile_g <- as.numeric(sum(d$qhat * alpha_g))

      var_quantile_g <- qlm_var_mean_bound(
        fhat = fhat_g,
        tau = d$tau,
        alpha = alpha_g,
        n_g = d$n
      )

      weights[[g]] <- alpha_g
      fhat[[g]] <- fhat_g

      if (d$type == "mean_quantile") {
        qlm_warn_mean_quantile_consistency(d, group = g)

        m_hat[g] <- d$mean
        var_m_hat[g] <- qlm_mean_variance_from_disclosure(
          d = d,
          default = var_quantile_g
        )
      } else {
        m_hat[g] <- m_quantile_g
        var_m_hat[g] <- var_quantile_g
      }
    }

    if (d$type == "mean") {
      m_hat[g] <- d$mean
      var_m_hat[g] <- qlm_mean_variance_from_disclosure(
        d = d,
        default = NA_real_
      )
      weights[[g]] <- NA_real_
      fhat[[g]] <- NA_real_
    }
  }

  if (any(!is.finite(m_hat))) {
    stop("Some group mean estimates are not finite.")
  }

  if (any(is.na(var_m_hat))) {
    warning(
      "Some groups do not have variance information for m_hat. ",
      "The coefficient covariance matrix and standard errors are set to NA."
    )
  }

  W <- diag(n_vec, nrow = G)
  A <- crossprod(X, W %*% X)

  qlm_warn_design_matrix(A)

  beta_hat <- solve(A, crossprod(X, W %*% m_hat))

  if (any(is.na(var_m_hat))) {
    vcov_beta <- matrix(NA_real_, nrow = p, ncol = p)
    se_beta <- rep(NA_real_, p)
  } else {
    meat <- matrix(0, nrow = p, ncol = p)

    for (g in seq_len(G)) {
      xg <- matrix(X[g, ], ncol = 1)
      meat <- meat + n_vec[g]^2 * var_m_hat[g] * (xg %*% t(xg))
    }

    A_inv <- solve(A)
    vcov_beta <- A_inv %*% meat %*% A_inv
    se_beta <- sqrt(diag(vcov_beta))
  }

  if (!is.null(colnames(X))) {
    names(beta_hat) <- colnames(X)
    names(se_beta) <- colnames(X)
    rownames(vcov_beta) <- colnames(X)
    colnames(vcov_beta) <- colnames(X)
  }

  out <- list(
    call = match.call(),
    coefficients = beta_hat,
    vcov = vcov_beta,
    se = se_beta,
    m_hat = m_hat,
    var_m_hat = var_m_hat,
    disclosure = disclosure,
    X = X,
    n_g = n_vec,
    weights = weights,
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
  cat("Quantile/mean disclosure linear model\n\n")
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
  cat("Quantile/mean disclosure linear model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}


as_qlm_disclosure <- function(disclosure = NULL,
                              qhat = NULL,
                              tau = NULL,
                              n_g = NULL) {
  if (!is.null(disclosure)) {
    if (is.data.frame(disclosure)) {
      disclosure <- qlm_disclosure_from_data_frame(disclosure)
    } else if (is.list(disclosure)) {
      if (length(disclosure) == 0) {
        stop("disclosure must be a non-empty list.")
      }
    } else {
      stop("disclosure must be either a data.frame or a list.")
    }

    for (g in seq_along(disclosure)) {
      disclosure[[g]] <- validate_qlm_disclosure_one(disclosure[[g]], group = g)
    }

    return(disclosure)
  }

  if (is.null(qhat) || is.null(tau) || is.null(n_g)) {
    stop("Either disclosure or qhat, tau, and n_g must be supplied.")
  }

  qhat <- as.matrix(qhat)
  tau <- as.numeric(tau)
  n_g <- as.numeric(n_g)

  if (!is.numeric(qhat)) {
    stop("qhat must be numeric.")
  }
  if (any(!is.finite(qhat))) {
    stop("qhat contains non-finite values.")
  }

  G <- nrow(qhat)
  K <- ncol(qhat)

  if (length(tau) != K) {
    stop("length(tau) must equal ncol(qhat).")
  }
  if (length(n_g) != G) {
    stop("length(n_g) must equal nrow(qhat).")
  }

  disclosure <- vector("list", G)

  for (g in seq_len(G)) {
    disclosure[[g]] <- list(
      type = "quantile",
      tau = tau,
      qhat = as.numeric(qhat[g, ]),
      n = n_g[g]
    )
    disclosure[[g]] <- validate_qlm_disclosure_one(disclosure[[g]], group = g)
  }

  disclosure
}


qlm_disclosure_from_data_frame <- function(df) {
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  if (nrow(df) == 0) {
    stop("disclosure data.frame must have at least one row.")
  }

  group_col <- qlm_find_first_existing_name(
    names(df),
    c("group", "group_id", "study", "study_id", "id")
  )

  if (is.null(group_col)) {
    stop(
      "disclosure data.frame must contain a group identifier column. ",
      "Accepted names are: group, group_id, study, study_id, or id."
    )
  }

  if (!"type" %in% names(df)) {
    stop("disclosure data.frame must contain a type column.")
  }

  if (!"n" %in% names(df)) {
    stop("disclosure data.frame must contain an n column.")
  }

  value_col <- qlm_find_first_existing_name(
    names(df),
    c("value", "qhat", "estimate")
  )

  has_mean_col <- "mean" %in% names(df)

  if (is.null(value_col) && !has_mean_col) {
    stop(
      "disclosure data.frame must contain either a value/qhat/estimate column ",
      "or a mean column."
    )
  }

  group_values <- unique(df[[group_col]])
  disclosure <- vector("list", length(group_values))

  for (j in seq_along(group_values)) {
    group_value <- group_values[j]
    dg <- df[df[[group_col]] == group_value, , drop = FALSE]

    type_vec <- as.character(dg$type)
    type_vec <- tolower(type_vec)

    is_quantile <- type_vec %in% c("quantile", "q", "quantiles")
    is_mean <- type_vec %in% c("mean", "m")

    if (!all(is_quantile | is_mean)) {
      stop(
        "Unsupported type found for group ", group_value, ". ",
        "Use 'quantile' or 'mean'."
      )
    }

    n_unique <- unique(stats::na.omit(as.numeric(dg$n)))

    if (length(n_unique) != 1) {
      stop("Each group must have exactly one non-missing value of n.")
    }

    n_g <- n_unique[1]

    out <- list(n = n_g)

    if (any(is_quantile)) {
      if (!"tau" %in% names(dg)) {
        stop("Quantile rows require a tau column.")
      }
      if (is.null(value_col)) {
        stop("Quantile rows require a value, qhat, or estimate column.")
      }

      q_rows <- dg[is_quantile, , drop = FALSE]
      tau_g <- as.numeric(q_rows$tau)
      qhat_g <- as.numeric(q_rows[[value_col]])

      ord <- order(tau_g)
      tau_g <- tau_g[ord]
      qhat_g <- qhat_g[ord]

      out$tau <- tau_g
      out$qhat <- qhat_g
    }

    if (any(is_mean)) {
      m_rows <- dg[is_mean, , drop = FALSE]

      if (has_mean_col && any(!is.na(m_rows$mean))) {
        mean_values <- unique(stats::na.omit(as.numeric(m_rows$mean)))
      } else if (!is.null(value_col)) {
        mean_values <- unique(stats::na.omit(as.numeric(m_rows[[value_col]])))
      } else {
        mean_values <- numeric(0)
      }

      if (length(mean_values) != 1) {
        stop("Each group with mean disclosure must have exactly one mean value.")
      }

      out$mean <- mean_values[1]

      if ("se" %in% names(m_rows)) {
        se_values <- unique(stats::na.omit(as.numeric(m_rows$se)))
        if (length(se_values) > 1) {
          stop("Each group must have at most one unique se value.")
        }
        if (length(se_values) == 1) {
          out$se <- se_values[1]
        }
      }

      if ("var" %in% names(m_rows)) {
        var_values <- unique(stats::na.omit(as.numeric(m_rows$var)))
        if (length(var_values) > 1) {
          stop("Each group must have at most one unique var value.")
        }
        if (length(var_values) == 1) {
          out$var <- var_values[1]
        }
      }

      if ("sd" %in% names(m_rows)) {
        sd_values <- unique(stats::na.omit(as.numeric(m_rows$sd)))
        if (length(sd_values) > 1) {
          stop("Each group must have at most one unique sd value.")
        }
        if (length(sd_values) == 1) {
          out$sd <- sd_values[1]
        }
      }
    }

    if (any(is_quantile) && any(is_mean)) {
      out$type <- "mean_quantile"
    } else if (any(is_quantile)) {
      out$type <- "quantile"
    } else {
      out$type <- "mean"
    }

    disclosure[[j]] <- out
  }

  names(disclosure) <- as.character(group_values)
  disclosure
}


qlm_find_first_existing_name <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) {
    return(NULL)
  }
  hit[1]
}


validate_qlm_disclosure_one <- function(d, group) {
  if (!is.list(d)) {
    stop("disclosure[[", group, "]] must be a list.")
  }

  if (is.null(d$type)) {
    stop("disclosure[[", group, "]] must contain type.")
  }

  d$type <- tolower(as.character(d$type))

  if (!d$type %in% c("quantile", "mean", "mean_quantile")) {
    stop(
      "disclosure[[", group, "]]$type must be one of ",
      "'quantile', 'mean', or 'mean_quantile'."
    )
  }

  if (is.null(d$n) ||
      length(d$n) != 1 ||
      !is.finite(d$n) ||
      d$n <= 0) {
    stop("disclosure[[", group, "]]$n must be a positive finite scalar.")
  }

  d$n <- as.numeric(d$n)

  if (d$type %in% c("quantile", "mean_quantile")) {
    if (is.null(d$tau) || is.null(d$qhat)) {
      stop(
        "disclosure[[", group, "]] with type '", d$type,
        "' must contain tau and qhat."
      )
    }

    d$tau <- as.numeric(d$tau)
    d$qhat <- as.numeric(d$qhat)

    if (length(d$tau) != length(d$qhat)) {
      stop(
        "In disclosure[[", group, "]], length(tau) must equal length(qhat)."
      )
    }
    if (length(d$tau) == 0) {
      stop("disclosure[[", group, "]]$tau must be non-empty.")
    }
    if (any(!is.finite(d$tau)) || any(!is.finite(d$qhat))) {
      stop(
        "tau and qhat must contain only finite values in disclosure[[",
        group, "]]."
      )
    }
    if (any(d$tau <= 0 | d$tau >= 1)) {
      stop("All tau values must lie in (0, 1) in disclosure[[", group, "]].")
    }
    if (is.unsorted(d$tau, strictly = TRUE)) {
      stop("tau must be strictly increasing in disclosure[[", group, "]].")
    }
    if (any(diff(d$qhat) < 0)) {
      stop("qhat must be weakly increasing in disclosure[[", group, "]].")
    }
  }

  if (d$type %in% c("mean", "mean_quantile")) {
    if (is.null(d$mean) ||
        length(d$mean) != 1 ||
        !is.finite(d$mean)) {
      stop(
        "Mean disclosure requires a finite scalar mean in disclosure[[",
        group, "]]."
      )
    }

    d$mean <- as.numeric(d$mean)

    if (!is.null(d$se)) {
      if (length(d$se) != 1 || !is.finite(d$se) || d$se < 0) {
        stop("disclosure[[", group, "]]$se must be a non-negative finite scalar.")
      }
      d$se <- as.numeric(d$se)
    }

    if (!is.null(d$var)) {
      if (length(d$var) != 1 || !is.finite(d$var) || d$var < 0) {
        stop("disclosure[[", group, "]]$var must be a non-negative finite scalar.")
      }
      d$var <- as.numeric(d$var)
    }

    if (!is.null(d$sd)) {
      if (length(d$sd) != 1 || !is.finite(d$sd) || d$sd < 0) {
        stop("disclosure[[", group, "]]$sd must be a non-negative finite scalar.")
      }
      d$sd <- as.numeric(d$sd)
    }
  }

  d
}


qlm_warn_disclosure_structure <- function(disclosure) {
  types <- vapply(disclosure, function(d) d$type, character(1))

  if (length(unique(types)) > 1) {
    message("Mixed disclosure detected: ", paste(unique(types), collapse = ", "), ".")
  }

  tau_lengths <- vapply(disclosure, function(d) {
    if (d$type %in% c("quantile", "mean_quantile")) {
      length(d$tau)
    } else {
      NA_integer_
    }
  }, integer(1))

  tau_lengths_nonmissing <- tau_lengths[!is.na(tau_lengths)]

  if (length(tau_lengths_nonmissing) > 1 &&
      length(unique(tau_lengths_nonmissing)) > 1) {
    warning(
      "Groups have different numbers of disclosed quantiles. ",
      "This is allowed, but precision may differ substantially across groups."
    )
  }

  for (g in seq_along(disclosure)) {
    d <- disclosure[[g]]

    if (d$type %in% c("quantile", "mean_quantile")) {
      max_gap <- max(diff(c(0, d$tau, 1)))

      if (max_gap > 0.5) {
        warning(
          "Group ", g, " has a coarse quantile grid with maximum tau gap ",
          round(max_gap, 3), ". Point inference may be sensitive to interpolation."
        )
      }

    }

    if (d$type == "mean" &&
        is.null(d$var) &&
        is.null(d$se) &&
        is.null(d$sd)) {
      warning(
        "Group ", g, " reports a mean but no se, var, or sd. ",
        "Its contribution to coefficient standard errors cannot be computed directly."
      )
    }
  }
}


qlm_warn_density <- function(fhat,
                             qhat,
                             tau,
                             group,
                             density_floor = 1e-8) {
  if (any(!is.finite(fhat))) {
    warning("Non-finite density estimate detected in group ", group, ".")
  }

  if (any(fhat >= 1 / density_floor)) {
    warning(
      "Density estimate in group ", group, " reached the numerical ceiling induced by ",
      "density_floor. This may indicate tied or nearly tied disclosed quantiles."
    )
  }

  dq <- diff(qhat)
  dtau <- diff(tau)

  if (any(dq <= density_floor * dtau)) {
    warning(
      "Group ", group, " has nearly tied adjacent disclosed quantiles. ",
      "Density and variance estimates may be unstable."
    )
  }
}


qlm_warn_mean_quantile_consistency <- function(d, group) {
  if (length(d$qhat) == 0) {
    return(invisible(NULL))
  }

  if (d$mean < min(d$qhat) || d$mean > max(d$qhat)) {
    warning(
      "In group ", group, ", the reported mean lies outside the range of disclosed quantiles. ",
      "This is not necessarily invalid, but it may indicate a heavy tail or inconsistent disclosure."
    )
  }

  invisible(NULL)
}


qlm_warn_design_matrix <- function(A) {
  if (qr(A)$rank < ncol(A)) {
    stop("The weighted design matrix is rank deficient.")
  }

  kappa_A <- kappa(A)

  if (is.finite(kappa_A) && kappa_A > 1e8) {
    warning(
      "The weighted design matrix appears ill-conditioned. ",
      "Coefficient estimates and standard errors may be numerically unstable."
    )
  }
}


qlm_mean_variance_from_disclosure <- function(d, default = NA_real_) {
  if (!is.null(d$var)) {
    return(as.numeric(d$var))
  }

  if (!is.null(d$se)) {
    return(as.numeric(d$se)^2)
  }

  if (!is.null(d$sd) && !is.null(d$n)) {
    return(as.numeric(d$sd)^2 / as.numeric(d$n))
  }

  default
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
      alpha[k] <- (tau[k] - tau[k - 1]) / 2 +
        (tau[k + 1] - tau[k]) / 2
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
