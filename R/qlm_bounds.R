#' Compute mean bounds from disclosed quantiles or observed means
#'
#' Compute sharp lower and upper bounds for group conditional means using
#' disclosed group-level quantiles and known support endpoints. If a group
#' reports an exact mean, the lower and upper bounds are both set equal to
#' the reported mean.
#'
#' @param qhat A G by K matrix of disclosed sample quantiles. Each row
#'   corresponds to a group and each column corresponds to a quantile level
#'   in `tau`. This is the legacy common-grid quantile-only interface.
#'   Ignored if `disclosure` is supplied.
#' @param tau A numeric vector of length K containing quantile levels in
#'   (0, 1), sorted in increasing order. This is the legacy common-grid
#'   quantile-only interface. Ignored if `disclosure` is supplied.
#' @param y_lower A numeric vector of length G, or a scalar, giving known
#'   lower support endpoints for each group. For quantile disclosures,
#'   support bounds must be supplied either through this argument and
#'   `y_upper`, or through columns `y_lower` and `y_upper` in `disclosure`.
#' @param y_upper A numeric vector of length G, or a scalar, giving known
#'   upper support endpoints for each group. For quantile disclosures,
#'   support bounds must be supplied either through `y_lower` and this
#'   argument, or through columns `y_lower` and `y_upper` in `disclosure`.
#' @param disclosure Optional group-level disclosure object. This may be
#'   either a list or a data.frame. If supplied, this takes precedence over
#'   `qhat` and `tau`.
#'
#'   A list input should contain one element per group. Each element should
#'   be a list with `type = "quantile"`, `type = "mean"`, or
#'   `type = "mean_quantile"`. Quantile disclosures require `tau` and
#'   `qhat`. Mean disclosures require `mean`. For quantile disclosures,
#'   support bounds must be supplied either through `y_lower` and `y_upper`
#'   inside each group element or through the function arguments `y_lower`
#'   and `y_upper`.
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
#'   For quantile disclosures, support bounds must be supplied either
#'   through columns `y_lower` and `y_upper` in the data.frame or through
#'   the function arguments `y_lower` and `y_upper`.
#'
#'#' @examples
#'
#' ## Legacy common-grid interface
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
#'   y_lower = 0,
#'   y_upper = 100
#' )
#'
#' bounds
#'
#' ## Disclosure data.frame interface
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
#'   \item{y_lower}{Lower support endpoints.}
#'   \item{y_upper}{Upper support endpoints.}
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
                       disclosure = NULL) {

  disclosure <- as_qlm_bounds_disclosure(
    disclosure = disclosure,
    qhat = qhat,
    tau = tau,
    y_lower = y_lower,
    y_upper = y_upper
  )

  G <- length(disclosure)

  m_lower <- numeric(G)
  m_upper <- numeric(G)
  y_lower_vec <- rep(NA_real_, G)
  y_upper_vec <- rep(NA_real_, G)
  weights_lower <- vector("list", G)
  weights_upper <- vector("list", G)

  for (g in seq_len(G)) {
    d <- disclosure[[g]]

    if (d$type %in% c("mean", "mean_quantile")) {
      m_lower[g] <- d$mean
      m_upper[g] <- d$mean
      y_lower_vec[g] <- if (!is.null(d$y_lower)) d$y_lower else NA_real_
      y_upper_vec[g] <- if (!is.null(d$y_upper)) d$y_upper else NA_real_
      weights_lower[[g]] <- NA_real_
      weights_upper[[g]] <- NA_real_
    }

    if (d$type == "quantile") {
      if (is.null(d$y_lower) || is.null(d$y_upper)) {
        stop(
          "Quantile-only group ", g,
          " requires y_lower and y_upper."
        )
      }

      bound_g <- qlm_bounds_one_group(
        qhat = d$qhat,
        tau = d$tau,
        y_lower = d$y_lower,
        y_upper = d$y_upper
      )

      m_lower[g] <- bound_g$m_lower
      m_upper[g] <- bound_g$m_upper
      y_lower_vec[g] <- d$y_lower
      y_upper_vec[g] <- d$y_upper
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
                                     y_upper = NULL) {
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
      y_upper = y_upper
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

  y_lower <- qlm_bounds_recycle_support(y_lower, G, "y_lower")
  y_upper <- qlm_bounds_recycle_support(y_upper, G, "y_upper")

  disclosure <- vector("list", G)

  for (g in seq_len(G)) {
    disclosure[[g]] <- list(
      type = "quantile",
      tau = tau,
      qhat = as.numeric(qhat[g, ]),
      y_lower = y_lower[g],
      y_upper = y_upper[g]
    )

    disclosure[[g]] <- validate_qlm_bounds_disclosure_one(
      disclosure[[g]],
      group = g
    )
  }

  disclosure
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

    if ("y_lower" %in% names(dg)) {
      y_lower_values <- unique(stats::na.omit(as.numeric(dg$y_lower)))
      if (length(y_lower_values) > 1) {
        stop("Each group must have at most one unique y_lower value.")
      }
      if (length(y_lower_values) == 1) {
        out$y_lower <- y_lower_values[1]
      }
    }

    if ("y_upper" %in% names(dg)) {
      y_upper_values <- unique(stats::na.omit(as.numeric(dg$y_upper)))
      if (length(y_upper_values) > 1) {
        stop("Each group must have at most one unique y_upper value.")
      }
      if (length(y_upper_values) == 1) {
        out$y_upper <- y_upper_values[1]
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
                                      y_upper = NULL) {
  G <- length(disclosure)

  if (!is.null(y_lower)) {
    y_lower <- qlm_bounds_recycle_support(y_lower, G, "y_lower")

    for (g in seq_len(G)) {
      if (is.null(disclosure[[g]]$y_lower)) {
        disclosure[[g]]$y_lower <- y_lower[g]
      }
    }
  }

  if (!is.null(y_upper)) {
    y_upper <- qlm_bounds_recycle_support(y_upper, G, "y_upper")

    for (g in seq_len(G)) {
      if (is.null(disclosure[[g]]$y_upper)) {
        disclosure[[g]]$y_upper <- y_upper[g]
      }
    }
  }

  disclosure
}


qlm_bounds_recycle_support <- function(x, G, name) {
  if (is.null(x)) {
    stop(name, " must be supplied for quantile-only mean bounds.")
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


validate_qlm_bounds_disclosure_one <- function(d, group) {
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
  }

  if (d$type == "quantile") {
    if (is.null(d$y_lower) || is.null(d$y_upper)) {
      stop(
        "Quantile-only disclosure requires y_lower and y_upper in disclosure[[",
        group, "]]."
      )
    }
  }

  if (!is.null(d$y_lower)) {
    if (length(d$y_lower) != 1 ||
        !is.finite(d$y_lower)) {
      stop("disclosure[[", group, "]]$y_lower must be a finite scalar.")
    }
    d$y_lower <- as.numeric(d$y_lower)
  }

  if (!is.null(d$y_upper)) {
    if (length(d$y_upper) != 1 ||
        !is.finite(d$y_upper)) {
      stop("disclosure[[", group, "]]$y_upper must be a finite scalar.")
    }
    d$y_upper <- as.numeric(d$y_upper)
  }

  if (!is.null(d$y_lower) && !is.null(d$y_upper)) {
    if (d$y_lower > d$y_upper) {
      stop(
        "In disclosure[[", group,
        "]], y_lower must be less than or equal to y_upper."
      )
    }
  }

  if (d$type == "quantile") {
    if (d$qhat[1] < d$y_lower) {
      stop(
        "In disclosure[[", group,
        "]], the first disclosed quantile must be greater than or equal to y_lower."
      )
    }

    if (d$qhat[length(d$qhat)] > d$y_upper) {
      stop(
        "In disclosure[[", group,
        "]], the last disclosed quantile must be less than or equal to y_upper."
      )
    }
  }

  d
}


qlm_bounds_one_group <- function(qhat, tau, y_lower, y_upper) {
  tau_full <- c(0, tau, 1)
  q_full <- c(y_lower, qhat, y_upper)

  delta_tau <- diff(tau_full)

  m_lower <- sum(delta_tau * q_full[-length(q_full)])
  m_upper <- sum(delta_tau * q_full[-1])

  weights_lower <- delta_tau[-length(delta_tau)]
  weights_upper <- delta_tau[-1]

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
#'   For quantile disclosures, support bounds must be supplied either
#'   through columns `y_lower` and `y_upper` in `disclosure` or through
#'   the function arguments `y_lower` and `y_upper`.
#' @param treatment Value of the `arm` column identifying treatment groups.
#'   For example, if `arm` takes values `"case"` and `"control"`,
#'   set `treatment = "case"`.
#' @param control Value of the `arm` column identifying control groups.
#'   For example, if `arm` takes values `"case"` and `"control"`,
#'   set `control = "control"`.
#' @param y_lower A scalar or numeric vector giving lower support endpoints.
#'   This argument is passed to `qlm_bounds()` for both treatment and control
#'   arms. For quantile disclosures, support bounds must be supplied either
#'   through this argument and `y_upper`, or through columns `y_lower` and
#'   `y_upper` in `disclosure`. If scalar, the same lower bound is used for
#'   all groups.
#' @param y_upper A scalar or numeric vector giving upper support endpoints.
#'   This argument is passed to `qlm_bounds()` for both treatment and control
#'   arms. For quantile disclosures, support bounds must be supplied either
#'   through `y_lower` and this argument, or through columns `y_lower` and
#'   `y_upper` in `disclosure`. If scalar, the same upper bound is used for
#'   all groups.
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
                                y_upper = NULL) {

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
    y_upper = y_upper
  )

  control_bounds <- qlm_bounds(
    disclosure = control_df,
    y_lower = y_lower,
    y_upper = y_upper
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
