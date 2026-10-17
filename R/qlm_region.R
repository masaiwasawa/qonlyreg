#' Simultaneous confidence region from quantile/mean disclosure
#'
#' Construct a simultaneous confidence region for regression coefficients
#' using disclosed group-level quantiles, observed means, or mixed disclosure.
#'
#' The function constructs estimated lower and upper bounds for each group
#' conditional mean, estimates the covariance matrix of these bounds, and
#' then forms a simultaneous confidence region for coefficients satisfying
#' \deqn{
#'   X_g'\beta \in [\underline m_g, \overline m_g].
#' }
#'
#' @param qhat A G by K matrix of disclosed sample quantiles. Each row
#'   corresponds to a group and each column corresponds to a quantile level
#'   in `tau`. This is the legacy common-grid quantile-only interface.
#'   Ignored if `disclosure` is supplied.
#' @param tau A numeric vector of length K containing quantile levels in
#'   (0, 1), sorted in increasing order. This is the legacy common-grid
#'   quantile-only interface. Ignored if `disclosure` is supplied.
#' @param X A G by p matrix of group-level covariates. Rows must correspond
#'   to the groups represented by `qhat` or `disclosure`, in the same order.
#' @param n_g A numeric vector of length G containing group sizes. This is
#'   the legacy common-grid quantile-only interface. Ignored if `disclosure`
#'   is supplied.
#' @param y_lower A numeric vector of length G, or a scalar, giving the
#'   exact lower endpoint of the conditional outcome support for each group,
#'   if known. When supplied, it is used both in constructing the group-mean
#'   bounds and as a deterministic boundary node in finite-difference
#'   estimation of quantile-function derivatives.
#' @param y_upper A numeric vector of length G, or a scalar, giving the
#'   exact upper endpoint of the conditional outcome support for each group,
#'   if known. When supplied, it is used both in constructing the group-mean
#'   bounds and as a deterministic boundary node in finite-difference
#'   estimation of quantile-function derivatives.
#' @param y_lower_bound A numeric vector of length G, or a scalar, giving a
#'   known lower bound on the conditional outcome support for each group.
#'   It is used to construct the group-mean bounds when `y_lower` is not
#'   supplied and is not used as a finite-difference boundary node.
#' @param y_upper_bound A numeric vector of length G, or a scalar, giving a
#'   known upper bound on the conditional outcome support for each group.
#'   It is used to construct the group-mean bounds when `y_upper` is not
#'   supplied and is not used as a finite-difference boundary node.
#' @param disclosure Optional group-level disclosure object. This may be
#'   either a list or a data.frame. If supplied, it takes precedence over
#'   `qhat`, `tau`, and `n_g`. Group-specific endpoint and support-bound
#'   values in `disclosure` take precedence over the corresponding
#'   function-level arguments.
#'
#'   A data.frame input should be in long format. It must contain a group
#'   identifier column, a `type` column, and an `n` column. The group
#'   identifier column is detected from one of `group`, `group_id`, `study`,
#'   `study_id`, or `id`. The `type` column indicates whether each row
#'   reports a quantile or a mean. Supported values are `quantile`, `q`,
#'   `quantiles`, `mean`, and `m`.
#'
#'   Quantile disclosures require `tau`, `qhat`, and `n`,
#'   together with a lower endpoint or support bound and an upper
#'   endpoint or support bound, supplied either within each group
#'   element or through the function arguments.
#'
#'   Mean rows require either a `mean` column or one disclosed-value column
#'   named `value`, `qhat`, or `estimate`. For inference with mean
#'   disclosures, at least one of `var`, `se`, or `sd` must also be supplied.
#'   Here `var` refers to the variance of the reported mean estimator, not
#'   the variance of the underlying observations.
#'
#'   For a quantile disclosure, each side of the support must be specified
#'   either by an exact endpoint (`y_lower` or `y_upper`) or by a support
#'   bound (`y_lower_bound` or `y_upper_bound`). These values may be supplied
#'   through columns in `disclosure` or through the corresponding function
#'   arguments.
#'
#'   A list input should contain one element per group. Each element should
#'   be a list with `type = "quantile"`, `type = "mean"`, or
#'   `type = "mean_quantile"`. Quantile disclosures require `tau`, `qhat`,
#'   `n`, `y_lower`, and `y_upper`, unless support bounds are supplied
#'   through the function arguments. Mean disclosures require `mean` and
#'   `n`, and at least one of `var`, `se`, or `sd`.
#' @param level Confidence level. Default is 0.95.
#' @param B Number of Gaussian simulation draws used to compute the
#'   simultaneous critical value.
#' @param density_floor A positive lower bound used to avoid numerical
#'   instability in density estimation.
#'
#' @details
#' For each group with a quantile disclosure, exact support endpoints may
#' be supplied through `y_lower` and `y_upper`, while support bounds that
#' need not be attained may be supplied through `y_lower_bound` and
#' `y_upper_bound`.
#'
#' Mean bounds use an exact endpoint when available and otherwise use the
#' corresponding support bound. Only exact endpoints supplied through
#' `y_lower` or `y_upper` are used as deterministic boundary nodes in
#' finite-difference estimation of quantile-function derivatives. Values
#' supplied through `y_lower_bound` and `y_upper_bound` do not enter the
#' finite-difference calculation.
#' The finite-difference derivative estimates are inverted to form the density
#' plug-ins used to construct the quantile-disclosure blocks of `Omega`.
#'
#' When a group reports a sample mean, the estimated lower and upper mean
#' bounds are both set equal to the reported mean. Its sampling uncertainty
#' is obtained from the supplied `var`, `se`, or `sd`.
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
#'   \item{const_mat}{Constraint matrix defining the simultaneous confidence
#'     region.}
#'   \item{const_dir}{Constraint directions.}
#'   \item{const_rhs}{Right-hand side vector of the constraints.}
#'   \item{disclosure}{Standardized group-level disclosure objects.}
#'   \item{X}{Group-level covariate matrix.}
#'   \item{n_g}{Group sizes.}
#'   \item{call}{Matched function call.}
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
#' region <- qlm_region(
#'   qhat = qhat,
#'   tau = tau,
#'   X = X,
#'   n_g = n_g,
#'   y_lower = 0,
#'   y_upper = 100,
#'   B = 100
#' )
#'
#' region
#' confint(region)
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
#' X <- cbind(
#'   intercept = 1,
#'   x = c(0, 1, 2)
#' )
#'
#' region <- qlm_region(
#'   disclosure = disclosure,
#'   X = X,
#'   y_lower = 0,
#'   y_upper = 100,
#'   B = 100
#' )
#'
#' region
#' confint(region)
#'
#' @export
qlm_region <- function(qhat = NULL,
                       tau = NULL,
                       X,
                       n_g = NULL,
                       y_lower = NULL,
                       y_upper = NULL,
                       disclosure = NULL,
                       level = 0.95,
                       B = 10000,
                       density_floor = 1e-8,
                       y_lower_bound = NULL,
                       y_upper_bound = NULL) {

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }
  if (any(!is.finite(X))) {
    stop("X contains non-finite values.")
  }
  if (!is.numeric(level) ||
      length(level) != 1 ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop("level must be a scalar in (0, 1).")
  }
  if (!is.numeric(B) ||
      length(B) != 1 ||
      !is.finite(B) ||
      B <= 0) {
    stop("B must be a positive scalar.")
  }
  if (!is.numeric(density_floor) ||
      length(density_floor) != 1 ||
      !is.finite(density_floor) ||
      density_floor <= 0) {
    stop("density_floor must be a positive finite scalar.")
  }

  B <- as.integer(B)

  disclosure <- as_qlm_region_disclosure(
    disclosure = disclosure,
    qhat = qhat,
    tau = tau,
    n_g = n_g,
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  G <- length(disclosure)
  N <- sum(vapply(disclosure, function(d) d$n, numeric(1)))

  if (nrow(X) != G) {
    stop("nrow(X) must equal the number of groups in disclosure.")
  }

  qlm_region_warn_disclosure_structure(disclosure)

  bounds <- qlm_region_bounds_from_disclosure(disclosure)

  Omega <- qlm_region_omega_from_disclosure(
    disclosure = disclosure,
    density_floor = density_floor
  )

  omega <- sqrt(diag(Omega))

  if (any(!is.finite(omega))) {
    stop("Estimated asymptotic standard deviations contain non-finite values.")
  }
  if (any(omega <= 0)) {
    stop(
      "Estimated asymptotic variances must be strictly positive. ",
      "This may occur when disclosed quantiles are nearly identical or ",
      "when mean disclosures have zero or missing variance."
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

  if (!is.null(colnames(X))) {
    colnames(const_mat) <- colnames(X)
  }

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
    disclosure = disclosure,
    X = X,
    n_g = vapply(disclosure, function(d) d$n, numeric(1)),
    call = match.call()
  )

  class(out) <- "qlm_region"
  out
}


#' @export
print.qlm_region <- function(x, ...) {
  cat("Quantile/mean disclosure simultaneous confidence region\n\n")
  cat("Confidence level:", x$level, "\n")
  cat("Critical value:", x$critical_value, "\n\n")
  cat("Linear constraints are stored in:\n")
  cat("  const_mat, const_dir, and const_rhs\n\n")
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
  cat("Quantile/mean disclosure simultaneous confidence region\n\n")

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
  cat("Quantile/mean disclosure coefficient-wise confidence intervals\n\n")
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
  cat("Quantile/mean disclosure coefficient-wise confidence intervals\n\n")

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
#' a simultaneous confidence region.
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

as_qlm_region_disclosure <- function(
    disclosure = NULL,
    qhat = NULL,
    tau = NULL,
    n_g = NULL,
    y_lower = NULL,
    y_upper = NULL,
    y_lower_bound = NULL,
    y_upper_bound = NULL) {

  # ---------------------------------------------------------------
  # 1. disclosure interface
  # ---------------------------------------------------------------
  if (!is.null(disclosure)) {
    if (is.data.frame(disclosure)) {
      disclosure <- qlm_region_disclosure_from_data_frame(disclosure)
    } else if (is.list(disclosure)) {
      if (length(disclosure) == 0L) {
        stop("disclosure must be a non-empty list.")
      }
    } else {
      stop("disclosure must be either a data.frame or a list.")
    }

    # Fill support information that is not already group-specific.
    disclosure <- qlm_region_attach_support(
      disclosure = disclosure,
      y_lower = y_lower,
      y_upper = y_upper,
      y_lower_bound = y_lower_bound,
      y_upper_bound = y_upper_bound
    )

    for (g in seq_along(disclosure)) {
      disclosure[[g]] <- validate_qlm_region_disclosure_one(
        disclosure[[g]],
        group = g
      )
    }

    return(disclosure)
  }

  # ---------------------------------------------------------------
  # 2. Legacy common-grid interface
  # ---------------------------------------------------------------
  if (is.null(qhat) ||
      is.null(tau) ||
      is.null(n_g)) {
    stop(
      "Either disclosure or qhat, tau, and n_g must be supplied."
    )
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

  # First construct the disclosure objects without requiring
  # either endpoints or bounds.
  disclosure <- vector("list", G)

  for (g in seq_len(G)) {
    disclosure[[g]] <- list(
      type = "quantile",
      tau = tau,
      qhat = as.numeric(qhat[g, ]),
      n = n_g[g]
    )
  }

  # Add whichever endpoint/bound arguments have been supplied.
  disclosure <- qlm_region_attach_support(
    disclosure = disclosure,
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  # The validator determines whether each group has sufficient
  # and mutually consistent support information.
  for (g in seq_along(disclosure)) {
    disclosure[[g]] <- validate_qlm_region_disclosure_one(
      disclosure[[g]],
      group = g
    )
  }

  disclosure
}


qlm_region_disclosure_from_data_frame <- function(df) {
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  if (nrow(df) == 0) {
    stop("disclosure data.frame must have at least one row.")
  }

  group_col <- qlm_region_find_first_existing_name(
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

  value_col <- qlm_region_find_first_existing_name(
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

    type_vec <- tolower(as.character(dg$type))

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

    out <- list(
      group_id = group_value,
      n = n_unique[1]
    )

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

    support_names <- c(
      "y_lower",
      "y_upper",
      "y_lower_bound",
      "y_upper_bound"
    )

    for (nm in support_names) {
      if (nm %in% names(dg)) {
        values <- unique(stats::na.omit(as.numeric(dg[[nm]])))
        if (length(values) > 1L) {
          stop("Each group must have at most one unique ", nm, " value.")
        }
        if (length(values) == 1L) {
          out[[nm]] <- values[1L]
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


qlm_region_attach_support <- function(
    disclosure,
    y_lower = NULL,
    y_upper = NULL,
    y_lower_bound = NULL,
    y_upper_bound = NULL) {

  G <- length(disclosure)

  support_values <- list(
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  for (nm in names(support_values)) {
    x <- support_values[[nm]]

    if (is.null(x)) {
      next
    }

    x <- qlm_region_recycle_vector(
      x = x,
      G = G,
      name = nm
    )

    for (g in seq_len(G)) {
      # Group-specific values in disclosure take precedence
      # over function-level arguments.
      if (is.null(disclosure[[g]][[nm]])) {
        disclosure[[g]][[nm]] <- x[g]
      }
    }
  }

  disclosure
}

validate_qlm_region_disclosure_one <- function(d, group) {
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

  if (is.null(d$n) || length(d$n) != 1L || !is.finite(d$n) || d$n <= 0) {
    stop(
      "disclosure[[", group, "]]$n must be a positive finite scalar."
    )
  }

  d$n <- as.numeric(d$n)

  # ---------------------------------------------------------------
  # Validate exact support endpoints and support bounds
  # ---------------------------------------------------------------
  support_names <- c(
    "y_lower",
    "y_upper",
    "y_lower_bound",
    "y_upper_bound"
  )

  for (nm in support_names) {
    if (!is.null(d[[nm]])) {
      if (length(d[[nm]]) != 1L) {
        stop(
          "disclosure[[", group, "]]$",
          nm,
          " must be a finite scalar."
        )
      }

      value <- suppressWarnings(as.numeric(d[[nm]]))

      if (!is.finite(value)) {
        stop(
          "disclosure[[", group, "]]$",
          nm,
          " must be a finite scalar."
        )
      }

      d[[nm]] <- value
    }
  }

  y_lower_endpoint <- d[["y_lower"]]
  y_upper_endpoint <- d[["y_upper"]]
  y_lower_bound <- d[["y_lower_bound"]]
  y_upper_bound <- d[["y_upper_bound"]]

  # A lower support bound cannot exceed the exact lower endpoint.
  if (!is.null(y_lower_bound) &&
      !is.null(y_lower_endpoint) &&
      y_lower_bound > y_lower_endpoint) {
    stop(
      "In disclosure[[", group,
      "]], y_lower_bound must not exceed y_lower."
    )
  }

  # An exact upper endpoint cannot exceed an upper support bound.
  if (!is.null(y_upper_endpoint) &&
      !is.null(y_upper_bound) &&
      y_upper_endpoint > y_upper_bound) {
    stop(
      "In disclosure[[", group,
      "]], y_upper must not exceed y_upper_bound."
    )
  }

  # Exact endpoints take precedence when constructing mean bounds.
  lower_for_bounds <- if (!is.null(y_lower_endpoint)) {
    y_lower_endpoint
  } else {
    y_lower_bound
  }

  upper_for_bounds <- if (!is.null(y_upper_endpoint)) {
    y_upper_endpoint
  } else {
    y_upper_bound
  }

  if (!is.null(lower_for_bounds) &&
      !is.null(upper_for_bounds) &&
      lower_for_bounds > upper_for_bounds) {
    stop(
      "In disclosure[[", group,
      "]], the maintained lower support value must not exceed ",
      "the maintained upper support value."
    )
  }

  # ---------------------------------------------------------------
  # Quantile disclosure
  # ---------------------------------------------------------------
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
        "In disclosure[[", group,
        "]], length(tau) must equal length(qhat)."
      )
    }

    if (length(d$tau) == 0L) {
      stop("disclosure[[", group, "]]$tau must be non-empty.")
    }

    if (any(!is.finite(d$tau)) ||
        any(!is.finite(d$qhat))) {
      stop(
        "tau and qhat must contain only finite values in disclosure[[",
        group, "]]."
      )
    }

    if (any(d$tau <= 0 | d$tau >= 1)) {
      stop(
        "All tau values must lie in (0, 1) in disclosure[[",
        group, "]]."
      )
    }

    if (is.unsorted(d$tau, strictly = TRUE)) {
      stop(
        "tau must be strictly increasing in disclosure[[",
        group, "]]."
      )
    }

    if (any(diff(d$qhat) < 0)) {
      stop(
        "qhat must be weakly increasing in disclosure[[",
        group, "]]."
      )
    }

    # Quantile-only disclosure requires information on both sides
    # to construct lower and upper mean bounds.
    if (d$type == "quantile" &&
        (is.null(lower_for_bounds) ||
         is.null(upper_for_bounds))) {
      stop(
        "Quantile disclosure requires a lower endpoint or lower ",
        "support bound and an upper endpoint or upper support bound ",
        "in disclosure[[", group, "]]."
      )
    }

    # Check compatibility of disclosed quantiles with the maintained
    # support interval. Exact endpoints are used when available;
    # otherwise the corresponding support bounds are used.
    if (!is.null(lower_for_bounds) &&
        d$qhat[1] < lower_for_bounds) {
      stop(
        "In disclosure[[", group,
        "]], the first disclosed quantile is below the maintained ",
        "lower support value."
      )
    }

    if (!is.null(upper_for_bounds) &&
        d$qhat[length(d$qhat)] > upper_for_bounds) {
      stop(
        "In disclosure[[", group,
        "]], the last disclosed quantile exceeds the maintained ",
        "upper support value."
      )
    }

    # With a single disclosed quantile, finite-difference covariance
    # estimation requires at least one exact support endpoint.
    # Support bounds alone are not finite-difference nodes.
    if (d$type == "quantile" &&
        length(d$tau) == 1L &&
        is.null(y_lower_endpoint) &&
        is.null(y_upper_endpoint)) {
      stop(
        "A quantile-only group with one disclosed quantile requires ",
        "at least one exact support endpoint for finite-difference ",
        "covariance estimation."
      )
    }
  }

  # ---------------------------------------------------------------
  # Mean disclosure
  # ---------------------------------------------------------------
  if (d$type %in% c("mean", "mean_quantile")) {
    if (is.null(d$mean) ||
        length(d$mean) != 1L ||
        !is.finite(d$mean)) {
      stop(
        "Mean disclosure requires a finite scalar mean in disclosure[[",
        group, "]]."
      )
    }

    d$mean <- as.numeric(d$mean)

    if (!is.null(d$se)) {
      if (length(d$se) != 1L ||
          !is.finite(d$se) ||
          d$se <= 0) {
        stop(
          "disclosure[[", group,
          "]]$se must be a positive finite scalar."
        )
      }
      d$se <- as.numeric(d$se)
    }

    if (!is.null(d$var)) {
      if (length(d$var) != 1L ||
          !is.finite(d$var) ||
          d$var <= 0) {
        stop(
          "disclosure[[", group,
          "]]$var must be a positive finite scalar."
        )
      }
      d$var <- as.numeric(d$var)
    }

    if (!is.null(d$sd)) {
      if (length(d$sd) != 1L ||
          !is.finite(d$sd) ||
          d$sd <= 0) {
        stop(
          "disclosure[[", group,
          "]]$sd must be a positive finite scalar."
        )
      }
      d$sd <- as.numeric(d$sd)
    }

    if (is.null(d$se) &&
        is.null(d$var) &&
        is.null(d$sd)) {
      stop(
        "Mean disclosure in qlm_region() requires uncertainty information: ",
        "provide se, var, or sd in disclosure[[", group, "]]."
      )
    }
  }

  d
}

qlm_region_lower_for_bounds <- function(d) {
  endpoint <- d[["y_lower"]]
  if (!is.null(endpoint)) { return(endpoint)}
  d[["y_lower_bound"]]
}

qlm_region_upper_for_bounds <- function(d) {
  endpoint <- d[["y_upper"]]
  if (!is.null(endpoint)) {return(endpoint)}
  d[["y_upper_bound"]]
}

qlm_region_bounds_from_disclosure <- function(disclosure) {
  G <- length(disclosure)

  m_lower <- numeric(G)
  m_upper <- numeric(G)

  # Exact endpoints, when supplied.
  y_lower_endpoint <- rep(NA_real_, G)
  y_upper_endpoint <- rep(NA_real_, G)

  # Values actually used to construct the mean bounds.
  # These may be exact endpoints or support bounds.
  y_lower_bound_used <- rep(NA_real_, G)
  y_upper_bound_used <- rep(NA_real_, G)

  weights_lower <- vector("list", G)
  weights_upper <- vector("list", G)

  for (g in seq_len(G)) {
    d <- disclosure[[g]]

    lower_for_bounds <- qlm_region_lower_for_bounds(d)
    upper_for_bounds <- qlm_region_upper_for_bounds(d)

    # Record exact endpoints separately.
    y_lower_endpoint_g <- d[["y_lower"]]
    y_upper_endpoint_g <- d[["y_upper"]]
    if (!is.null(y_lower_endpoint_g)) {y_lower_endpoint[g] <- y_lower_endpoint_g}
    if (!is.null(y_upper_endpoint_g)) {y_upper_endpoint[g] <- y_upper_endpoint_g}

    # Record the values actually used for mean bounding.
    if (!is.null(lower_for_bounds)) {
      y_lower_bound_used[g] <- lower_for_bounds
    }

    if (!is.null(upper_for_bounds)) {
      y_upper_bound_used[g] <- upper_for_bounds
    }

    if (d$type == "quantile") {
      # This should already have been guaranteed by
      # validate_qlm_region_disclosure_one().
      if (is.null(lower_for_bounds) ||
          is.null(upper_for_bounds)) {
        stop(
          "Internal error: quantile disclosure group ",
          g,
          " does not have support information on both sides."
        )
      }

      bound_g <- qlm_region_bounds_one_quantile_group(
        qhat = d$qhat,
        tau = d$tau,
        y_lower_bound = lower_for_bounds,
        y_upper_bound = upper_for_bounds
      )

      m_lower[g] <- bound_g$m_lower
      m_upper[g] <- bound_g$m_upper
      weights_lower[[g]] <- bound_g$weights_lower
      weights_upper[[g]] <- bound_g$weights_upper
    }

    if (d$type %in% c("mean", "mean_quantile")) {
      m_lower[g] <- d$mean
      m_upper[g] <- d$mean
      weights_lower[[g]] <- NA_real_
      weights_upper[[g]] <- NA_real_
    }
  }

  out <- list(
    m_lower = m_lower,
    m_upper = m_upper,

    # Exact endpoints only; NA when only a support bound was supplied.
    y_lower = y_lower_endpoint,
    y_upper = y_upper_endpoint,

    # Effective values used to construct the mean bounds.
    y_lower_bound = y_lower_bound_used,
    y_upper_bound = y_upper_bound_used,

    weights_lower = weights_lower,
    weights_upper = weights_upper
  )

  class(out) <- "qlm_bounds"
  out
}

qlm_region_bounds_one_quantile_group <- function(
    qhat,
    tau,
    y_lower_bound,
    y_upper_bound) {

  tau_full <- c(0, tau, 1)
  q_full <- c(y_lower_bound, qhat, y_upper_bound)

  delta_tau <- diff(tau_full)

  m_lower <- sum(
    delta_tau * q_full[-length(q_full)]
  )

  m_upper <- sum(
    delta_tau * q_full[-1]
  )

  # Coefficients on the disclosed sample quantiles.
  weights_lower <- delta_tau[-1]
  weights_upper <- delta_tau[-length(delta_tau)]

  list(
    m_lower = m_lower,
    m_upper = m_upper,
    weights_lower = weights_lower,
    weights_upper = weights_upper
  )
}


qlm_region_omega_from_disclosure <- function(disclosure,
                                             density_floor = 1e-8) {
  G <- length(disclosure)
  n_vec <- vapply(disclosure, function(d) d$n, numeric(1))
  N <- sum(n_vec)

  Omega <- matrix(0, nrow = 2 * G, ncol = 2 * G)

  for (g in seq_len(G)) {
    d <- disclosure[[g]]

    if (d$type == "quantile") {
      Omega_g <- qlm_region_omega_one_quantile_group(
        qhat = d$qhat,
        tau = d$tau,
        n_g = d$n,
        N = N,
        y_lower = d[["y_lower"]],
        y_upper = d[["y_upper"]],
        density_floor = density_floor
      )
    }

    if (d$type %in% c("mean", "mean_quantile")) {
      var_mean_g <- qlm_region_mean_variance_from_disclosure(d)
      Omega_g <- N * var_mean_g * matrix(1, nrow = 2, ncol = 2)
    }

    idx <- c(2 * g - 1, 2 * g)
    Omega[idx, idx] <- Omega_g
  }

  Omega
}


qlm_region_omega_one_quantile_group <- function(qhat,
                                                tau,
                                                n_g,
                                                N,
                                                y_lower,
                                                y_upper,
                                                density_floor = 1e-8) {
  tau_full <- c(0, tau, 1)
  delta_tau <- diff(tau_full)

  aL <- delta_tau[-1]
  aU <- delta_tau[-length(delta_tau)]

  fhat <- qlm_density_at_quantiles_with_endpoints(
    qhat = qhat,
    tau = tau,
    y_lower = y_lower,
    y_upper = y_upper,
    density_floor = density_floor
  )

  cov_tau <- outer(tau, tau, function(a, b) pmin(a, b) - a * b)
  denom <- outer(fhat, fhat, "*")
  Sigma_g <- cov_tau / denom

  (N / n_g) *
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
}


qlm_density_at_quantiles_with_endpoints <- function(
    qhat,
    tau,
    y_lower = NULL,
    y_upper = NULL,
    density_floor = 1e-8) {

  qhat <- as.numeric(qhat)
  tau <- as.numeric(tau)

  K <- length(tau)
  qprime_hat <- numeric(K)

  if (K == 1L) {
    if (!is.null(y_lower) && !is.null(y_upper)) {
      qprime_hat[1] <- y_upper - y_lower
    } else if (!is.null(y_lower)) {
      qprime_hat[1] <-
        (qhat[1] - y_lower) / tau[1]
    } else if (!is.null(y_upper)) {
      qprime_hat[1] <-
        (y_upper - qhat[1]) / (1 - tau[1])
    } else {
      stop(
        "One disclosed quantile requires at least one exact support ",
        "endpoint for finite-difference density estimation."
      )
    }
  } else {
    if (!is.null(y_lower)) {
      qprime_hat[1] <-
        (qhat[2] - y_lower) / tau[2]
    } else {
      qprime_hat[1] <-
        (qhat[2] - qhat[1]) /
        (tau[2] - tau[1])
    }

    if (K >= 3L) {
      for (k in 2:(K - 1L)) {
        qprime_hat[k] <-
          (qhat[k + 1L] - qhat[k - 1L]) /
          (tau[k + 1L] - tau[k - 1L])
      }
    }

    if (!is.null(y_upper)) {
      qprime_hat[K] <-
        (y_upper - qhat[K - 1L]) /
        (1 - tau[K - 1L])
    } else {
      qprime_hat[K] <-
        (qhat[K] - qhat[K - 1L]) /
        (tau[K] - tau[K - 1L])
    }
  }

  qprime_hat <- pmax(qprime_hat, density_floor)
  fhat <- 1 / qprime_hat

  pmax(fhat, density_floor)
}


qlm_region_mean_variance_from_disclosure <- function(d) {
  if (!is.null(d$var)) {
    return(as.numeric(d$var))
  }

  if (!is.null(d$se)) {
    return(as.numeric(d$se)^2)
  }

  if (!is.null(d$sd) && !is.null(d$n)) {
    return(as.numeric(d$sd)^2 / as.numeric(d$n))
  }

  stop("Mean disclosure requires se, var, or sd.")
}


qlm_region_warn_disclosure_structure <- function(disclosure) {
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

    if (d$type %in% c("quantile")) {
      max_gap <- max(diff(c(0, d$tau, 1)))

      if (max_gap > 0.5) {
        warning(
          "Group ", g, " has a coarse quantile grid with maximum tau gap ",
          round(max_gap, 3), ". The confidence region may be wide."
        )
      }

      lower_endpoint <- d[["y_lower"]]
      upper_endpoint <- d[["y_upper"]]

      if (length(d$tau) == 1L) {
        if (!is.null(lower_endpoint) && !is.null(upper_endpoint)) {
          endpoint_description <- "both exact support endpoints"
        } else if (!is.null(lower_endpoint)) {
          endpoint_description <- "the exact lower support endpoint"
        } else if (!is.null(upper_endpoint)) {
          endpoint_description <- "the exact upper support endpoint"
        } else {
          stop(
            "Internal error: a quantile-only group with one disclosed ",
            "quantile reached covariance estimation without an exact ",
            "support endpoint."
          )
        }

        warning(
          "Group ", g, " has only one disclosed quantile. ",
          "The finite-difference density plug-in uses ",
          endpoint_description,
          "."
        )
      }
    }
  }
}


qlm_region_recycle_vector <- function(x, G, name) {
  if (is.null(x)) {
    stop(name, " must be supplied.")
  }

  if (length(x) == 1) {
    x <- rep(x, G)
  }

  x <- as.numeric(x)

  if (length(x) != G) {
    stop(name, " must be a scalar or a vector of length G.")
  }

  if (any(!is.finite(x))) {
    stop(name, " contains non-finite values.")
  }

  x
}


qlm_region_find_first_existing_name <- function(x, candidates) {
  hit <- candidates[candidates %in% x]

  if (length(hit) == 0) {
    return(NULL)
  }

  hit[1]
}





#' Simultaneous confidence region for treatment-control contrasts
#'
#' Construct a simultaneous confidence region for regression coefficients
#' when the outcome is a treatment-control mean contrast. This function is
#' intended for two-arm studies, meta-analysis applications, and related
#' settings where treatment and control summaries are disclosed separately.
#'
#' For each matched group `g`, the contrast is
#' \deqn{
#'   \theta_g = m_{g,T} - m_{g,C}.
#' }
#' The function constructs bounds for each \eqn{\theta_g}, estimates their
#' covariance matrix, and then forms a simultaneous confidence region for
#' coefficients satisfying
#' \deqn{
#'   X_g'\beta \in [\underline\theta_g, \overline\theta_g].
#' }
#'
#' @param disclosure A data.frame containing arm-level disclosure data.
#'   The data.frame must be in long format and contain:
#'   \describe{
#'     \item{group identifier column}{A column identifying matched studies,
#'       pairs, experiments, or other units. The column name is detected
#'       automatically from one of `group`, `group_id`, `study`,
#'       `study_id`, or `id`.}
#'     \item{arm}{A column named `arm` identifying treatment and control
#'       arms.}
#'     \item{type}{A column indicating whether each row reports a quantile
#'       or a mean. Supported values are `quantile`, `q`, `quantiles`,
#'       `mean`, and `m`.}
#'     \item{n}{A column giving the arm-level sample size.}
#'   }
#'
#'   Quantile rows require a `tau` column and one disclosed-value column.
#'   The disclosed-value column may be named `value`, `qhat`, or `estimate`.
#'
#'   Mean rows require either a `mean` column or one disclosed-value column
#'   named `value`, `qhat`, or `estimate`. For inference with mean
#'   disclosures, at least one of `var`, `se`, or `sd` must also be supplied.
#'   Here `var` refers to the variance of the reported mean estimator, not
#'   the variance of the underlying observations.
#'
#'   For quantile disclosures, each side must be specified by either an
#'   exact endpoint (`y_lower` or `y_upper`) or a support bound
#'   (`y_lower_bound` or `y_upper_bound`). These values may be supplied
#'   through columns in `disclosure` or through the corresponding
#'   function arguments.
#'
#' @param X A G by p matrix of group-level covariates. Rows must correspond
#'   to the matched groups after alignment by the detected group identifier.
#' @param treatment Value of the `arm` column identifying treatment groups.
#'   For example, if `arm` takes values `"case"` and `"control"`, set
#'   `treatment = "case"`.
#' @param control Value of the `arm` column identifying control groups.
#'   For example, if `arm` takes values `"case"` and `"control"`, set
#'   `control = "control"`.
#' @param y_lower A scalar or numeric vector giving exact lower support
#'   endpoints, if known. Exact endpoints are used both in mean bounding
#'   and as finite-difference boundary nodes.
#' @param y_upper A scalar or numeric vector giving exact upper support
#'   endpoints, if known. Exact endpoints are used both in mean bounding
#'   and as finite-difference boundary nodes.
#' @param y_lower_bound A scalar or numeric vector giving known lower
#'   support bounds. These values are used in mean bounding when
#'   `y_lower` is not supplied and are not used as finite-difference nodes.
#' @param y_upper_bound A scalar or numeric vector giving known upper
#'   support bounds. These values are used in mean bounding when
#'   `y_upper` is not supplied and are not used as finite-difference nodes.
#' @param level Confidence level. Default is 0.95.
#' @param B Number of Gaussian simulation draws used to compute the
#'   simultaneous critical value.
#' @param density_floor A positive lower bound used to avoid numerical
#'   instability in density estimation.
#'
#' @return
#' An object of class `"qlm_contrast_region"` and `"qlm_region"` containing:
#' \describe{
#'   \item{bounds}{Estimated contrast bounds.}
#'   \item{Omega}{Estimated covariance matrix of the estimated contrast bounds.}
#'   \item{critical_value}{Simulated simultaneous critical value.}
#'   \item{level}{Confidence level.}
#'   \item{lower_adjusted}{Adjusted lower constraints.}
#'   \item{upper_adjusted}{Adjusted upper constraints.}
#'   \item{const_mat}{Constraint matrix defining the simultaneous confidence region.}
#'   \item{const_dir}{Constraint directions.}
#'   \item{const_rhs}{Right-hand side vector of the constraints.}
#'   \item{id}{Matched group identifiers.}
#'   \item{group_col}{Detected group identifier column.}
#'   \item{treatment}{Treatment arm value.}
#'   \item{control}{Control arm value.}
#'   \item{X}{Group-level covariate matrix.}
#'   \item{n_g}{Total sample size within each matched treatment-control pair.}
#'   \item{call}{Matched function call.}
#' }
#'
#' @examples
#' disclosure <- data.frame(
#'   study = c(
#'     1, 1, 1,
#'     1, 1, 1,
#'     2,
#'     2
#'   ),
#'   arm = c(
#'     "case", "case", "case",
#'     "control", "control", "control",
#'     "case",
#'     "control"
#'   ),
#'   type = c(
#'     "quantile", "quantile", "quantile",
#'     "quantile", "quantile", "quantile",
#'     "mean",
#'     "mean"
#'   ),
#'   tau = c(
#'     0.25, 0.50, 0.75,
#'     0.25, 0.50, 0.75,
#'     NA,
#'     NA
#'   ),
#'   value = c(
#'     10, 15, 20,
#'     8, 12, 16,
#'     25,
#'     20
#'   ),
#'   n = c(
#'     100, 100, 100,
#'     100, 100, 100,
#'     120,
#'     120
#'   ),
#'   se = c(
#'     NA, NA, NA,
#'     NA, NA, NA,
#'     1.5,
#'     1.2
#'   )
#' )
#'
#' X <- matrix(1, nrow = 2, ncol = 1)
#' colnames(X) <- "intercept"
#'
#' region <- qlm_contrast_region(
#'   disclosure = disclosure,
#'   X = X,
#'   treatment = "case",
#'   control = "control",
#'   y_lower = 0,
#'   y_upper = 100,
#'   B = 100
#' )
#'
#' confint(region)
#'
#' @export
qlm_contrast_region <- function(disclosure,
                                X,
                                treatment,
                                control,
                                y_lower = NULL,
                                y_upper = NULL,
                                level = 0.95,
                                B = 10000,
                                density_floor = 1e-8,
                                y_lower_bound = NULL,
                                y_upper_bound = NULL) {

  if (!is.data.frame(disclosure)) {
    stop("disclosure must be a data.frame.")
  }
  if (nrow(disclosure) == 0) {
    stop("disclosure must have at least one row.")
  }

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }
  if (any(!is.finite(X))) {
    stop("X contains non-finite values.")
  }
  if (!is.numeric(level) || length(level) != 1 ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be a scalar in (0, 1).")
  }
  if (!is.numeric(B) || length(B) != 1 ||
      !is.finite(B) || B <= 0) {
    stop("B must be a positive scalar.")
  }
  if (!is.numeric(density_floor) || length(density_floor) != 1 ||
      !is.finite(density_floor) || density_floor <= 0) {
    stop("density_floor must be a positive finite scalar.")
  }

  B <- as.integer(B)

  group_col <- qlm_region_find_first_existing_name(
    names(disclosure),
    c("group", "group_id", "study", "study_id", "id")
  )

  if (is.null(group_col)) {
    stop(
      "disclosure must contain a group identifier column. ",
      "Accepted names are: group, group_id, study, study_id, or id."
    )
  }

  if (!"arm" %in% names(disclosure)) {
    stop("disclosure must contain an arm column.")
  }

  if (missing(treatment)) {
    stop("treatment must be supplied.")
  }
  if (missing(control)) {
    stop("control must be supplied.")
  }
  if (identical(treatment, control)) {
    stop("treatment and control must be different arm values.")
  }

  treat_df <- disclosure[disclosure$arm == treatment, , drop = FALSE]
  control_df <- disclosure[disclosure$arm == control, , drop = FALSE]

  if (nrow(treat_df) == 0) {
    stop("No treatment-arm rows found.")
  }
  if (nrow(control_df) == 0) {
    stop("No control-arm rows found.")
  }

  treat_ids <- unique(as.character(treat_df[[group_col]]))
  control_ids <- unique(as.character(control_df[[group_col]]))

  missing_control <- setdiff(treat_ids, control_ids)
  missing_treat <- setdiff(control_ids, treat_ids)

  if (length(missing_control) > 0) {
    stop(
      "Some treatment ids have no matching control ids: ",
      paste(missing_control, collapse = ", ")
    )
  }

  if (length(missing_treat) > 0) {
    stop(
      "Some control ids have no matching treatment ids: ",
      paste(missing_treat, collapse = ", ")
    )
  }

  matched_ids <- sort(unique(treat_ids))

  treat_df <- treat_df[
    order(match(as.character(treat_df[[group_col]]), matched_ids),
          seq_len(nrow(treat_df))),
    ,
    drop = FALSE
  ]

  control_df <- control_df[
    order(match(as.character(control_df[[group_col]]), matched_ids),
          seq_len(nrow(control_df))),
    ,
    drop = FALSE
  ]

  treat_disclosure <- qlm_region_disclosure_from_data_frame(treat_df)
  control_disclosure <- qlm_region_disclosure_from_data_frame(control_df)

  treat_disclosure <- qlm_region_attach_support(
    disclosure = treat_disclosure,
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  control_disclosure <- qlm_region_attach_support(
    disclosure = control_disclosure,
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  for (g in seq_along(treat_disclosure)) {
    treat_disclosure[[g]] <- validate_qlm_region_disclosure_one(
      treat_disclosure[[g]],
      group = g
    )
  }

  for (g in seq_along(control_disclosure)) {
    control_disclosure[[g]] <- validate_qlm_region_disclosure_one(
      control_disclosure[[g]],
      group = g
    )
  }

  G <- length(treat_disclosure)

  if (length(control_disclosure) != G) {
    stop("Treatment and control disclosures must have the same number of groups.")
  }

  if (nrow(X) != G) {
    stop("nrow(X) must equal the number of matched treatment-control groups.")
  }

  n_treat <- vapply(treat_disclosure, function(d) d$n, numeric(1))
  n_control <- vapply(control_disclosure, function(d) d$n, numeric(1))
  n_pair <- n_treat + n_control
  N <- sum(n_pair)

  bounds_treat <- qlm_region_bounds_from_disclosure(treat_disclosure)
  bounds_control <- qlm_region_bounds_from_disclosure(control_disclosure)

  Omega_treat <- qlm_region_omega_from_disclosure(
    disclosure = treat_disclosure,
    density_floor = density_floor
  )

  Omega_control <- qlm_region_omega_from_disclosure(
    disclosure = control_disclosure,
    density_floor = density_floor
  )

  N_treat <- sum(n_treat)
  N_control <- sum(n_control)

  Omega_treat <- Omega_treat * (N / N_treat)
  Omega_control <- Omega_control * (N / N_control)

  theta_lower <- bounds_treat$m_lower - bounds_control$m_upper
  theta_upper <- bounds_treat$m_upper - bounds_control$m_lower

  if (any(theta_lower > theta_upper)) {
    stop("Some lower contrast bounds exceed upper contrast bounds.")
  }

  Omega <- matrix(0, nrow = 2 * G, ncol = 2 * G)

  for (g in seq_len(G)) {
    idx <- c(2 * g - 1, 2 * g)

    Omega_T_g <- Omega_treat[idx, idx, drop = FALSE]
    Omega_C_g <- Omega_control[idx, idx, drop = FALSE]

    Omega_g <- matrix(NA_real_, nrow = 2, ncol = 2)

    ## theta_lower = treatment_lower - control_upper
    ## theta_upper = treatment_upper - control_lower
    Omega_g[1, 1] <- Omega_T_g[1, 1] + Omega_C_g[2, 2]
    Omega_g[2, 2] <- Omega_T_g[2, 2] + Omega_C_g[1, 1]
    Omega_g[1, 2] <- Omega_T_g[1, 2] + Omega_C_g[2, 1]
    Omega_g[2, 1] <- Omega_g[1, 2]

    Omega[idx, idx] <- Omega_g
  }

  omega <- sqrt(diag(Omega))

  if (any(!is.finite(omega))) {
    stop("Estimated asymptotic standard deviations contain non-finite values.")
  }
  if (any(omega <= 0)) {
    stop(
      "Estimated asymptotic variances must be strictly positive. ",
      "This may occur when disclosed quantiles are nearly identical or ",
      "when mean disclosures have zero or missing variance."
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

  lower_adj <- theta_lower - adj[seq(1, 2 * G, by = 2)]
  upper_adj <- theta_upper + adj[seq(2, 2 * G, by = 2)]

  const_mat <- rbind(X, X)
  const_dir <- c(rep(">=", G), rep("<=", G))
  const_rhs <- c(lower_adj, upper_adj)

  if (!is.null(colnames(X))) {
    colnames(const_mat) <- colnames(X)
  }

  bounds <- list(
    theta_lower = theta_lower,
    theta_upper = theta_upper
  )
  class(bounds) <- "qlm_contrast_bounds"

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
    id = matched_ids,
    group_col = group_col,
    treatment = treatment,
    control = control,
    treatment_disclosure = treat_disclosure,
    control_disclosure = control_disclosure,
    X = X,
    n_g = n_pair,
    call = match.call()
  )

  class(out) <- c("qlm_contrast_region", "qlm_region")
  out
}


#' @export
print.qlm_contrast_region <- function(x, ...) {
  cat("Treatment-control contrast simultaneous confidence region\n\n")
  cat("Confidence level:", x$level, "\n")
  cat("Critical value:", x$critical_value, "\n\n")
  cat("Linear constraints are stored in:\n")
  cat("  const_mat, const_dir, and const_rhs\n\n")
  cat("Number of constraints:", nrow(x$const_mat), "\n")
  cat("Number of coefficients:", ncol(x$const_mat), "\n")
  invisible(x)
}
