#' Compute DUH and Related Heterogeneity Measures
#'
#' @description
#' Computes the Dominance-Unevenness-Heterogeneity (DUH) index along with the
#' Gini-Simpson Index (GSI) and Shannon Entropy (SE) for each observation.
#' The inner evenness calculation uses high-precision arithmetic via \pkg{Rmpfr}.
#'
#' @param data       A data frame.
#' @param groups     A length-1 string prefix identifying group proportion columns,
#'                   or a character vector of exact column names.
#' @param evepar     Numeric. Evenness parameter; must be strictly > 1. Default 2.
#' @param round_digits Integer. Decimal places for output columns. Default 4.
#' @param show_groups  Logical. Print identified group names? Default TRUE.
#' @param prefix     Character or NULL. Prefix prepended to DUH, GSI, SE output columns.
#' @param precision  Integer. Bit-precision passed to \code{Rmpfr::mpfr()}. Default 128.
#'
#' @return The original data frame with columns appended:
#'   \code{sigma_1}, \code{evenness}, \code{DUH} (or \code{<prefix>DUH}),
#'   \code{GSI} (or \code{<prefix>GSI}), \code{SE} (or \code{<prefix>SE}).
#'
#' @details
#' Group proportion columns must sum to 1 per row. If they do not, automatic
#' row-wise normalisation is applied with a warning.
#'
#' High-precision arithmetic (\pkg{Rmpfr}) is used only for the \code{sigma_tilde}
#' and evenness (\code{psi}) calculation to avoid catastrophic cancellation when
#' group proportions are near-equal.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(g_a = c(0.5, 0.8), g_b = c(0.3, 0.1), g_c = c(0.2, 0.1))
#' get_duh(df, groups = "g_", evepar = 2)
#' }
#'
#' @importFrom tidyverse tidyverse_packages
#' @importFrom Rmpfr mpfr
#' @export
get_duh <- function(data,
                    groups,
                    evepar      = 2,
                    round_digits = 4,
                    show_groups  = TRUE,
                    prefix       = NULL,
                    precision    = 128) {

  if (!requireNamespace("tidyverse", quietly = TRUE))
    stop("Package 'tidyverse' is required. Install with: install.packages('tidyverse')")
  if (!requireNamespace("Rmpfr", quietly = TRUE))
    stop("Package 'Rmpfr' is required for high-precision math. Install with: install.packages('Rmpfr')")

  library(tidyverse)
  library(Rmpfr)

  if (evepar <= 1)
    stop("evepar must be strictly greater than 1 for convexity")

  # ── Identify group columns ──────────────────────────────────────────────────
  if (length(groups) > 1) {
    g2      <- groups
    abs_s   <- length(groups)
    subdata <- data[, groups, drop = FALSE]
  } else {
    col_names <- colnames(data)
    g         <- stringr::str_detect(col_names, paste0("^", groups))
    abs_s     <- sum(g)
    g2        <- col_names[g]
    subdata   <- data[, g, drop = FALSE]
  }

  if (abs_s < 2)
    stop("At least 2 group columns are required.")

  # ── Check / normalise proportions ───────────────────────────────────────────
  row_sums <- rowSums(subdata)
  if (!isTRUE(all.equal(round(row_sums, 0), rep(1, nrow(subdata))))) {
    message("Total Group Proportions Do Not Add Up To 1, automatic scaling applied")
    subdata <- subdata / row_sums
  } else {
    message("Check: Group Proportions Add Up to 1")
  }

  # ── High-precision evenness (psi) for a single row ─────────────────────────
  # Uses Rmpfr to avoid catastrophic cancellation in the sigma_tilde calculation.
  compute_psi <- function(minor_p_vec) {
    if (sum(minor_p_vec) == 0) return(0)           # edge: all weight on sigma_1
    mp_evepar    <- mpfr(evepar,       precBits = precision)
    mp_factor    <- log(mp_evepar)
    mp_minor     <- mpfr(minor_p_vec,  precBits = precision)
    mp_minor_sum <- sum(mp_minor)
    sigma_tilde  <- mp_factor * (mp_minor / mp_minor_sum)
    psi_terms    <- abs(sigma_tilde - mp_factor / (abs_s - 1))^mp_evepar
    as.numeric(1 - sum(psi_terms)^(1 / mp_evepar) / mp_factor)
  }

  # ── Row-wise calculation ────────────────────────────────────────────────────
  n            <- nrow(subdata)
  sigma_1_vec  <- numeric(n)
  evenness_vec <- numeric(n)
  duh_vec      <- numeric(n)
  gsi_vec      <- numeric(n)
  se_vec       <- numeric(n)

  for (i in seq_len(n)) {
    p_vals   <- as.numeric(subdata[i, ])
    sorted_p <- sort(p_vals, decreasing = TRUE)
    sigma_1  <- sorted_p[1]
    minor_p  <- sorted_p[-1]

    # Evenness (high precision)
    psi <- compute_psi(minor_p)

    # DUH: log(sigma_1) = 0 when sigma_1 = 1 (monopoly) → DUH = 0
    duh <- if (sigma_1 <= 0 || sigma_1 >= 1) {
      0
    } else {
      val <- -(log(sigma_1) / log(abs_s)) * psi
      ifelse(is.nan(val), 0, val)
    }

    # HHI / GSI
    hhi <- sum(p_vals^2)
    gsi <- 1 - hhi

    # Shannon Entropy
    plp <- ifelse(p_vals > 0, -p_vals * log(p_vals), 0)
    se  <- sum(plp)

    sigma_1_vec[i]  <- sigma_1
    evenness_vec[i] <- round(psi, round_digits)
    duh_vec[i]      <- round(duh, round_digits)
    gsi_vec[i]      <- round(gsi, round_digits)
    se_vec[i]       <- round(se,  round_digits)
  }

  # ── Assemble output ─────────────────────────────────────────────────────────
  duh_col <- if (!is.null(prefix)) paste0(prefix, "DUH") else "DUH"
  gsi_col <- if (!is.null(prefix)) paste0(prefix, "GSI") else "GSI"
  se_col  <- if (!is.null(prefix)) paste0(prefix, "SE")  else "SE"

  data$sigma_1    <- sigma_1_vec
  data$evenness   <- evenness_vec
  data[[duh_col]] <- duh_vec
  data[[gsi_col]] <- gsi_vec
  data[[se_col]]  <- se_vec

  message(sprintf("Universe has %d categories", abs_s))
  if (show_groups) message(paste(g2, collapse = ", "))
  message(sprintf("Using %s-metric", evepar))

  return(data)
}
