#' Compute mean bounds from disclosed quantiles or observed means
#'
#' Compute sharp lower and upper bounds for group conditional means.
#' For quantile-only groups, the bounds use disclosed group-level
#' quantiles and known lower and upper bounds on the conditional
#' outcome support. If a group reports an exact mean, the lower and
#' upper bounds are both set equal to the reported mean.
#'
#' @param qhat A G by K matrix of disclosed sample quantiles. Each row
#'   corresponds to a group and each column corresponds to a quantile level
#'   in `tau`. This is the legacy common-grid quantile-only interface.
#'   Ignored if `disclosure` is supplied.
#' @param tau A numeric vector of length K containing quantile levels in
#'   (0, 1), sorted in increasing order. This is the legacy common-grid
#'   quantile-only interface. Ignored if `disclosure` is supplied.
#' @param y_lower A numeric vector of length G, or a scalar, giving the
#'   exact lower endpoint of the conditional outcome support for each group,
#'   if known. When supplied, it is used as the lower support bound in mean
#'   bounding.
#' @param y_upper A numeric vector of length G, or a scalar, giving the
#'   exact upper endpoint of the conditional outcome support for each group,
#'   if known. When supplied, it is used as the upper support bound in mean
#'   bounding.
#' @param y_lower_bound A numeric vector of length G, or a scalar, giving a
#'   known lower bound on the conditional outcome support. This argument is
#'   used only when `y_lower` is not supplied.
#' @param y_upper_bound A numeric vector of length G, or a scalar, giving a
#'   known upper bound on the conditional outcome support. This argument is
#'   used only when `y_upper` is not supplied.
#' @param disclosure Optional group-level disclosure object. This may be
#'   either a list or a data.frame. If supplied, this takes precedence over
#'   `qhat` and `tau`.
#'
#'   A list input should contain one element per group. Each element should
#'   be a list with `type = "quantile"`, `type = "mean"`, or
#'   `type = "mean_quantile"`. Quantile disclosures require `tau` and
#'   `qhat`. Mean disclosures require `mean`.
#'
#'   For a quantile-only group, the lower side of the support must be
#'   specified through either `y_lower`, for an exact lower support
#'   endpoint, or `y_lower_bound`, for a lower support bound. Similarly,
#'   the upper side must be specified through either `y_upper` or
#'   `y_upper_bound`. These values may be included in each group element
#'   or supplied through the corresponding function arguments.
#'
#'   A data.frame input should be in long format. It must contain a group
#'   identifier column and a `type` column. The group identifier column is
#'   detected from one of `group`, `group_id`, `study`, `study_id`, or `id`.
#'   The `type` column indicates whether each row reports a quantile or a
#'   mean. Supported values are `quantile`, `q`, `quantiles`, `mean`, and
#'   `m`.
#'
#'   Quantile rows require a `tau` column and one disclosed-value column.
#'   The disclosed-value column may be named `value`, `qhat`, or `estimate`.
#'
#'   Mean rows require either a `mean` column or one disclosed-value column
#'   named `value`, `qhat`, or `estimate`.
#'
#'   For quantile-only groups, support information may be supplied through
#'   columns `y_lower`, `y_upper`, `y_lower_bound`, and
#'   `y_upper_bound`, or through the corresponding function arguments.
#'   For each group, either `y_lower` or `y_lower_bound` is required on
#'   the lower side, and either `y_upper` or `y_upper_bound` is required
#'   on the upper side.
#'
#' @details
#' For a quantile-only group, mean bounds require a known value
#' on each side of the outcome support. An exact endpoint may be supplied
#' through `y_lower` or `y_upper`; otherwise a support bound may be supplied
#' through `y_lower_bound` or `y_upper_bound`.
#'
#' When both an endpoint and a bound are supplied on the same side, the
#' endpoint is used because it is the sharper restriction. The supplied
#' values must be mutually consistent.
#'
#' For a group with an exact mean disclosure, the reported mean determines
#' both the lower and upper mean bounds, so support bounds are not needed
#' to bound that group's mean.
#'
#' @examples
#'
#' ## Legacy common-grid interface using support bounds
#' qhat <- rbind(
#'   c(10, 15, 20),
#'   c(20, 25, 30)
#' )
#'
#' tau <- c(0.25, 0.50, 0.75)
#'
#' bounds <- qlm_bounds(
#'   qhat = qhat,
#'   tau = tau,
#'   y_lower_bound = 0,
#'   y_upper_bound = 100
#' )
#'
#' bounds
#'
#' ## Disclosure data.frame interface using exact support endpoints
#' disclosure <- data.frame(
#'   group = c(
#'     1, 1, 1,
#'     2
#'   ),
#'   type = c(
#'     "quantile", "quantile", "quantile",
#'     "mean"
#'   ),
#'   tau = c(
#'     0.25, 0.50, 0.75,
#'     NA
#'   ),
#'   value = c(
#'     10, 15, 20,
#'     25
#'   )
#' )
#'
#' bounds <- qlm_bounds(
#'   disclosure = disclosure,
#'   y_lower = 0,
#'   y_upper = 100
#' )
#'
#' bounds
#'
#' @return
#' An object of class `"qlm_bounds"` containing:
#' \describe{
#'   \item{m_lower}{Estimated lower bounds for group conditional means.}
#'   \item{m_upper}{Estimated upper bounds for group conditional means.}
#'   \item{disclosure}{Standardized group-level disclosure objects.}
#'   \item{qhat}{Input matrix of disclosed sample quantiles for the legacy
#'     common-grid interface. Set to `NULL` for heterogeneous disclosure.}
#'   \item{tau}{Input quantile levels for the legacy common-grid interface.
#'     Set to `NULL` for heterogeneous disclosure.}
#'   \item{y_lower}{Exact lower support endpoints, where supplied.}
#'   \item{y_upper}{Exact upper support endpoints, where supplied.}
#'   \item{y_lower_bound}{Lower support bounds, where supplied.}
#'   \item{y_upper_bound}{Upper support bounds, where supplied.}
#'   \item{weights_lower}{A list of weights on disclosed quantiles for lower
#'     bounds.}
#'   \item{weights_upper}{A list of weights on disclosed quantiles for upper
#'     bounds.}
#' }
#'
#' @export
qlm_bounds <- function(qhat = NULL,
                       tau = NULL,
                       y_lower = NULL,
                       y_upper = NULL,
                       disclosure = NULL,
                       y_lower_bound = NULL,
                       y_upper_bound = NULL) {

  disclosure <- as_qlm_bounds_disclosure(
    disclosure = disclosure,
    qhat = qhat,
    tau = tau,
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  G <- length(disclosure)

  m_lower <- numeric(G)
  m_upper <- numeric(G)
  y_lower_vec <- rep(NA_real_, G)
  y_upper_vec <- rep(NA_real_, G)
  y_lower_bound_vec <- rep(NA_real_, G)
  y_upper_bound_vec <- rep(NA_real_, G)
  weights_lower <- vector("list", G)
  weights_upper <- vector("list", G)

  for (g in seq_len(G)) {
    d <- disclosure[[g]]

    if (d$type %in% c("mean", "mean_quantile")) {
      m_lower[g] <- d$mean
      m_upper[g] <- d$mean
      y_lower_endpoint <- d[["y_lower"]]
      y_upper_endpoint <- d[["y_upper"]]
      y_lower_vec[g] <- if (!is.null(y_lower_endpoint)) {y_lower_endpoint} else { NA_real_}
      y_upper_vec[g] <- if (!is.null(y_upper_endpoint)) {y_upper_endpoint} else { NA_real_}
      y_lower_bound_vec[g] <- if (!is.null(d$y_lower_bound)) d$y_lower_bound else NA_real_
      y_upper_bound_vec[g] <- if (!is.null(d$y_upper_bound)) d$y_upper_bound else NA_real_
      weights_lower[[g]] <- NA_real_
      weights_upper[[g]] <- NA_real_
    }

    if (d$type == "quantile") {
      lower_for_bounds <- qlm_bounds_lower_for_bounds(d)
      upper_for_bounds <- qlm_bounds_upper_for_bounds(d)

      if (is.null(lower_for_bounds) || is.null(upper_for_bounds)) {
        stop(
          "Quantile-only group ", g,
          " requires a lower endpoint or lower support bound and ",
          "an upper endpoint or upper support bound."
        )
      }

      bound_g <- qlm_bounds_one_group(
        qhat = d$qhat,
        tau = d$tau,
        y_lower_bound = lower_for_bounds,
        y_upper_bound = upper_for_bounds
      )

      m_lower[g] <- bound_g$m_lower
      m_upper[g] <- bound_g$m_upper
      y_lower_endpoint <- d[["y_lower"]]
      y_upper_endpoint <- d[["y_upper"]]
      y_lower_vec[g] <- if (!is.null(y_lower_endpoint)) {y_lower_endpoint} else { NA_real_}
      y_upper_vec[g] <- if (!is.null(y_upper_endpoint)) {y_upper_endpoint} else { NA_real_}
      y_lower_bound_vec[g] <- if (!is.null(d$y_lower_bound)) d$y_lower_bound else NA_real_
      y_upper_bound_vec[g] <- if (!is.null(d$y_upper_bound)) d$y_upper_bound else NA_real_
      weights_lower[[g]] <- bound_g$weights_lower
      weights_upper[[g]] <- bound_g$weights_upper
    }
  }

  if (any(m_lower > m_upper)) {
    stop("Some lower mean bounds exceed upper mean bounds.")
  }

  legacy_qhat <- NULL
  legacy_tau <- NULL

  if (is.null(disclosure[[1]]$group_id) && !is.null(qhat) && !is.null(tau)) {
    legacy_qhat <- as.matrix(qhat)
    legacy_tau <- as.numeric(tau)
  }

  out <- list(
    m_lower = m_lower,
    m_upper = m_upper,
    disclosure = disclosure,
    qhat = legacy_qhat,
    tau = legacy_tau,
    y_lower = y_lower_vec,
    y_upper = y_upper_vec,
    y_lower_bound = y_lower_bound_vec,
    y_upper_bound = y_upper_bound_vec,
    weights_lower = weights_lower,
    weights_upper = weights_upper
  )

  class(out) <- "qlm_bounds"
  out
}


#' @export
print.qlm_bounds <- function(x, ...) {
  cat("Quantile/mean disclosure mean bounds\n\n")
  out <- data.frame(
    m_lower = x$m_lower,
    m_upper = x$m_upper
  )
  print(out)
  invisible(x)
}


as_qlm_bounds_disclosure <- function(disclosure = NULL,
                                     qhat = NULL,
                                     tau = NULL,
                                     y_lower = NULL,
                                     y_upper = NULL,
                                     y_lower_bound = NULL,
                                     y_upper_bound = NULL) {
  if (!is.null(disclosure)) {
    if (is.data.frame(disclosure)) {
      disclosure <- qlm_bounds_disclosure_from_data_frame(disclosure)
    } else if (is.list(disclosure)) {
      if (length(disclosure) == 0) {
        stop("disclosure must be a non-empty list.")
      }
    } else {
      stop("disclosure must be either a data.frame or a list.")
    }

    disclosure <- qlm_bounds_attach_support(
      disclosure = disclosure,
      y_lower = y_lower,
      y_upper = y_upper,
      y_lower_bound = y_lower_bound,
      y_upper_bound = y_upper_bound
    )

    for (g in seq_along(disclosure)) {
      disclosure[[g]] <- validate_qlm_bounds_disclosure_one(
        disclosure[[g]],
        group = g
      )
    }

    return(disclosure)
  }

  if (is.null(qhat) || is.null(tau)) {
    stop("Either disclosure or qhat and tau must be supplied.")
  }

  qhat <- as.matrix(qhat)
  tau <- as.numeric(tau)

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

  y_lower <- qlm_bounds_recycle_optional(y_lower, G, "y_lower")
  y_upper <- qlm_bounds_recycle_optional(y_upper, G, "y_upper")
  y_lower_bound <- qlm_bounds_recycle_optional(
    y_lower_bound, G, "y_lower_bound"
  )
  y_upper_bound <- qlm_bounds_recycle_optional(
    y_upper_bound, G, "y_upper_bound"
  )

  disclosure <- vector("list", G)

  for (g in seq_len(G)) {
    d <- list(
      type = "quantile",
      tau = tau,
      qhat = as.numeric(qhat[g, ])
    )

    if (!is.null(y_lower)) {
      d$y_lower <- y_lower[g]
    }
    if (!is.null(y_upper)) {
      d$y_upper <- y_upper[g]
    }
    if (!is.null(y_lower_bound)) {
      d$y_lower_bound <- y_lower_bound[g]
    }
    if (!is.null(y_upper_bound)) {
      d$y_upper_bound <- y_upper_bound[g]
    }

    disclosure[[g]] <- validate_qlm_bounds_disclosure_one(
      d,
      group = g
    )
  }

  disclosure
}


qlm_bounds_lower_for_bounds <- function(d) {
  if (!is.null(d[["y_lower"]])) {return(d[["y_lower"]])}
  d$y_lower_bound
}

qlm_bounds_upper_for_bounds <- function(d) {
  if (!is.null(d[["y_upper"]])) {return(d[["y_upper"]])}
  d$y_upper_bound
}

qlm_bounds_disclosure_from_data_frame <- function(df) {
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  if (nrow(df) == 0) {
    stop("disclosure data.frame must have at least one row.")
  }

  group_col <- qlm_bounds_find_first_existing_name(
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

  value_col <- qlm_bounds_find_first_existing_name(
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

    out <- list(group_id = group_value)

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
    }

    support_names <- c(
      "y_lower",
      "y_upper",
      "y_lower_bound",
      "y_upper_bound"
    )

    for (nm in support_names) {
      if (nm %in% names(dg)) {
        values <- unique(
          stats::na.omit(
            as.numeric(dg[[nm]])
          )
        )

        if (length(values) > 1L) {
          stop(
            "Each group must have at most one unique ",
            nm,
            " value."
          )
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


qlm_bounds_attach_support <- function(disclosure,
                                      y_lower = NULL,
                                      y_upper = NULL,
                                      y_lower_bound = NULL,
                                      y_upper_bound = NULL) {

  G <- length(disclosure)

  supplied <- list(
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  for (nm in names(supplied)) {
    x <- supplied[[nm]]

    if (is.null(x)) {next}
    x <- qlm_bounds_recycle_optional(x, G, nm)

    for (g in seq_len(G)) {
      if (is.null(disclosure[[g]][[nm]])) {
        disclosure[[g]][[nm]] <- x[g]
      }
    }
  }

  disclosure
}

qlm_bounds_recycle_optional <- function(x, G, name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (length(x) == 1L) {
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

validate_qlm_bounds_disclosure_one <- function(d, group) {

  # ---------------------------------------------------------------
  # 1. Basic disclosure structure
  # ---------------------------------------------------------------

  if (!is.list(d)) {
    stop("disclosure[[", group, "]] must be a list.")
  }

  if (is.null(d$type) || length(d$type) != 1L) {
    stop(
      "disclosure[[", group,
      "]] must contain a scalar type."
    )
  }

  d$type <- tolower(as.character(d$type))

  if (is.na(d$type) ||
      !(d$type %in% c("quantile", "mean", "mean_quantile"))) {
    stop(
      "disclosure[[", group, "]]$type must be one of ",
      "'quantile', 'mean', or 'mean_quantile'."
    )
  }

  # ---------------------------------------------------------------
  # 2. Quantile disclosure
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
      stop(
        "disclosure[[", group,
        "]]$tau must be non-empty."
      )
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
  }

  # ---------------------------------------------------------------
  # 3. Mean disclosure
  # ---------------------------------------------------------------

  if (d$type %in% c("mean", "mean_quantile")) {

    if (is.null(d$mean) || length(d$mean) != 1L) {
      stop(
        "Mean disclosure requires a finite scalar mean in disclosure[[",
        group, "]]."
      )
    }

    d$mean <- suppressWarnings(as.numeric(d$mean))

    if (!is.finite(d$mean)) {
      stop(
        "Mean disclosure requires a finite scalar mean in disclosure[[",
        group, "]]."
      )
    }
  }

  # ---------------------------------------------------------------
  # 4. Support endpoints and support bounds
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
          "disclosure[[", group, "]]$", nm,
          " must be a finite scalar."
        )
      }

      value <- suppressWarnings(as.numeric(d[[nm]]))

      if (!is.finite(value)) {
        stop(
          "disclosure[[", group, "]]$", nm,
          " must be a finite scalar."
        )
      }

      d[[nm]] <- value
    }
  }

  # A lower support bound cannot exceed the exact lower endpoint.
  if (!is.null(d[["y_lower"]]) &&
      !is.null(d$y_lower_bound) &&
      d$y_lower_bound > d[["y_lower"]]) {
    stop(
      "In disclosure[[", group,
      "]], y_lower_bound must not exceed y_lower."
    )
  }

  # An exact upper endpoint cannot exceed an upper support bound.
  if (!is.null(d[["y_upper"]]) &&
      !is.null(d$y_upper_bound) &&
      d[["y_upper"]] > d$y_upper_bound) {
    stop(
      "In disclosure[[", group,
      "]], y_upper must not exceed y_upper_bound."
    )
  }

  lower_for_bounds <- qlm_bounds_lower_for_bounds(d)
  upper_for_bounds <- qlm_bounds_upper_for_bounds(d)

  # Quantile-only disclosure requires enough information to bound
  # both tails. Mean and mean-quantile disclosures do not, because
  # their mean is directly reported.
  if (d$type == "quantile") {

    if (is.null(lower_for_bounds) ||
        is.null(upper_for_bounds)) {
      stop(
        "Quantile-only disclosure requires a lower endpoint or lower ",
        "support bound and an upper endpoint or upper support bound in ",
        "disclosure[[", group, "]]."
      )
    }
  }

  if (!is.null(lower_for_bounds) &&
      !is.null(upper_for_bounds) &&
      lower_for_bounds > upper_for_bounds) {
    stop(
      "In disclosure[[", group,
      "]], the lower support value must not exceed the upper support value."
    )
  }

  # ---------------------------------------------------------------
  # 5. Compatibility of quantiles with support information
  # ---------------------------------------------------------------

  if (d$type %in% c("quantile", "mean_quantile")) {

    if (!is.null(lower_for_bounds) &&
        d$qhat[1L] < lower_for_bounds) {
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
  }

  # ---------------------------------------------------------------
  # 6. Compatibility of a reported mean with support information
  # ---------------------------------------------------------------

  if (d$type %in% c("mean", "mean_quantile")) {

    if (!is.null(lower_for_bounds) &&
        d$mean < lower_for_bounds) {
      stop(
        "In disclosure[[", group,
        "]], the reported mean is below the maintained lower ",
        "support value."
      )
    }

    if (!is.null(upper_for_bounds) &&
        d$mean > upper_for_bounds) {
      stop(
        "In disclosure[[", group,
        "]], the reported mean exceeds the maintained upper ",
        "support value."
      )
    }
  }

  d
}


qlm_bounds_one_group <- function(
    qhat,
    tau,
    y_lower_bound,
    y_upper_bound) {

  tau_full <- c(0, tau, 1)
  q_full <- c(y_lower_bound, qhat, y_upper_bound)

  delta_tau <- diff(tau_full)

  m_lower <- sum(delta_tau * q_full[-length(q_full)])
  m_upper <- sum(delta_tau * q_full[-1])

  weights_lower <- delta_tau[-1L]
  weights_upper <- delta_tau[-length(delta_tau)]

  list(
    m_lower = m_lower,
    m_upper = m_upper,
    weights_lower = weights_lower,
    weights_upper = weights_upper
  )
}


qlm_bounds_find_first_existing_name <- function(x, candidates) {
  hit <- candidates[candidates %in% x]

  if (length(hit) == 0) {
    return(NULL)
  }

  hit[1]
}


#' Compute contrast bounds from treatment and control disclosures
#'
#' Compute lower and upper bounds for treatment-control mean contrasts
#' using arm-level disclosure data. This function is intended for settings
#' such as two-arm studies or meta-analysis applications where each study
#' reports summary information separately for treatment and control groups.
#'
#' Internally, the function first splits the disclosure data into treatment
#' and control arms, applies `qlm_bounds()` separately to each arm, aligns
#' the resulting bounds by group identifier, and then computes sharp bounds
#' for the treatment-control contrast.
#'
#' For each matched group `g`, the contrast is defined as
#' \deqn{
#'   \theta_g = m_{g,T} - m_{g,C}.
#' }
#' Given bounds
#' \deqn{
#'   m_{g,T} \in [\underline m_{g,T}, \overline m_{g,T}]
#' }
#' and
#' \deqn{
#'   m_{g,C} \in [\underline m_{g,C}, \overline m_{g,C}],
#' }
#' the resulting contrast bounds are
#' \deqn{
#'   \underline\theta_g
#'   =
#'   \underline m_{g,T} - \overline m_{g,C}
#' }
#' and
#' \deqn{
#'   \overline\theta_g
#'   =
#'   \overline m_{g,T} - \underline m_{g,C}.
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
#'       `mean`, and `m`, following the same convention as `qlm_bounds()`.}
#'   }
#'
#'   Quantile rows require a `tau` column and one disclosed-value column.
#'   The disclosed-value column may be named `value`, `qhat`, or `estimate`.
#'
#'   Mean rows require either a `mean` column or one disclosed-value column
#'   named `value`, `qhat`, or `estimate`.
#'
#'   For quantile-only disclosures, the lower side of the support must be
#'   specified through either `y_lower` or `y_lower_bound`, and the upper
#'   side through either `y_upper` or `y_upper_bound`. These values may be
#'   supplied through columns in `disclosure` or through the corresponding
#'   function arguments.
#' @param treatment Value of the `arm` column identifying treatment groups.
#'   For example, if `arm` takes values `"case"` and `"control"`,
#'   set `treatment = "case"`.
#' @param control Value of the `arm` column identifying control groups.
#'   For example, if `arm` takes values `"case"` and `"control"`,
#'   set `control = "control"`.
#' @param y_lower A scalar or numeric vector giving exact lower support
#'   endpoints, if known. This argument is passed to `qlm_bounds()` for
#'   both treatment and control arms. If scalar, the same lower endpoint
#'   is used for all groups.
#' @param y_upper A scalar or numeric vector giving exact upper support
#'   endpoints, if known. This argument is passed to `qlm_bounds()` for
#'   both treatment and control arms. If scalar, the same upper endpoint
#'   is used for all groups.
#' @param y_lower_bound A scalar or numeric vector giving known lower
#'   bounds on the conditional outcome support when the exact lower
#'   endpoints are not known. This argument is passed to `qlm_bounds()`
#'   for both treatment and control arms. If scalar, the same lower bound
#'   is used for all groups.
#' @param y_upper_bound A scalar or numeric vector giving known upper
#'   bounds on the conditional outcome support when the exact upper
#'   endpoints are not known. This argument is passed to `qlm_bounds()`
#'   for both treatment and control arms. If scalar, the same upper bound
#'   is used for all groups.
#'
#' @return
#' An object of class `"qlm_contrast_bounds"` containing:
#' \describe{
#'   \item{theta_lower}{Estimated lower bounds for treatment-control mean
#'     contrasts.}
#'   \item{theta_upper}{Estimated upper bounds for treatment-control mean
#'     contrasts.}
#'   \item{id}{Matched group identifiers used to align treatment and control
#'     arms.}
#'   \item{group_col}{Name of the group identifier column detected from
#'     `disclosure`.}
#'   \item{treatment}{The treatment arm value used in the function call.}
#'   \item{control}{The control arm value used in the function call.}
#'   \item{treatment_bounds}{The treatment-arm `"qlm_bounds"` object.}
#'   \item{control_bounds}{The control-arm `"qlm_bounds"` object.}
#'   \item{call}{The matched function call.}
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
#'   )
#' )
#'
#' contrast <- qlm_contrast_bounds(
#'   disclosure = disclosure,
#'   treatment = "case",
#'   control = "control",
#'   y_lower = 0,
#'   y_upper = 100
#' )
#'
#' contrast
#'
#' @export
qlm_contrast_bounds <- function(disclosure,
                                treatment,
                                control,
                                y_lower = NULL,
                                y_upper = NULL,
                                y_lower_bound = NULL,
                                y_upper_bound = NULL) {

  if (!is.data.frame(disclosure)) {
    stop("disclosure must be a data.frame.")
  }

  if (nrow(disclosure) == 0) {
    stop("disclosure must have at least one row.")
  }

  group_col <- qlm_bounds_find_first_existing_name(
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

  treat_ids <- unique(treat_df[[group_col]])
  control_ids <- unique(control_df[[group_col]])

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
    order(match(treat_df[[group_col]], matched_ids)),
    ,
    drop = FALSE
  ]

  control_df <- control_df[
    order(match(control_df[[group_col]], matched_ids)),
    ,
    drop = FALSE
  ]

  treat_bounds <- qlm_bounds(
    disclosure = treat_df,
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  control_bounds <- qlm_bounds(
    disclosure = control_df,
    y_lower = y_lower,
    y_upper = y_upper,
    y_lower_bound = y_lower_bound,
    y_upper_bound = y_upper_bound
  )

  stopifnot(
    length(treat_bounds$m_lower) ==
      length(control_bounds$m_lower)
  )

  theta_lower <- treat_bounds$m_lower - control_bounds$m_upper
  theta_upper <- treat_bounds$m_upper - control_bounds$m_lower

  if (any(theta_lower > theta_upper)) {
    stop("Some lower contrast bounds exceed upper contrast bounds.")
  }

  out <- list(
    theta_lower = theta_lower,
    theta_upper = theta_upper,
    id = matched_ids,
    group_col = group_col,
    treatment = treatment,
    control = control,
    treatment_bounds = treat_bounds,
    control_bounds = control_bounds,
    call = match.call()
  )

  class(out) <- "qlm_contrast_bounds"
  out
}


#' @export
print.qlm_contrast_bounds <- function(x, ...) {
  cat("Treatment-control contrast bounds\n\n")

  out <- data.frame(
    id = x$id,
    theta_lower = x$theta_lower,
    theta_upper = x$theta_upper
  )

  print(out)
  invisible(x)
}
