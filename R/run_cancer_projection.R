# ==========================================================
# Project expected cancer cases for a cohort (FULL SCRIPT)
# - strict age-band parsing & validation (catches "o" vs "0")
# - overlap & coverage checks (single-age level)
# - supports 85+ via open-ended band extended to age_max
# - decimals kept internally; rounding only for display
# ==========================================================

#' Project expected cancer cases for a cohort
#'
#' @description
#' Estimates expected cancer cases over time for a defined cohort, using flat
#' age-band incidence and mortality rates. Expands 5-year age-band counts to
#' single-year ages and simulates annual aging, incidence, and mortality.
#'
#' Separate tables of age-band rates must be supplied for cancer incidence and
#' all-cause mortality, for both females and males. Tables include `band`
#' and `rate_per100k`. Open-ended bands like "85+" are supported and extended
#' to `age_max`.
#'
#' @param counts_5y Tibble/data.frame with columns `age_band` (e.g., "40-44")
#'   and `N` (combined sexes baseline counts).
#' @param female_share,male_share Numeric values between 0 and 1; internally normalized
#'   if they do not sum to 1.
#' @param female_inc_wide,female_mort_5y Tibbles with `band`, `rate_per100k`
#'   for female incidence and mortality.
#' @param male_inc_wide,male_mort_5y Same for male.
#' @param study_start,study_end Integers, recruitment window.
#' @param age_min,age_max Integers defining modeled age range.
#' @param begin_year Optional integer; default = floor(midpoint of
#'   `study_start` and `study_end`).
#' @param end_year Final projection year (inclusive).
#' @param diag_years Optional integer vector of years to return in diagnostics;
#'   default includes `begin_year`, 2030, 2040 (kept within range).
#' @param print_markdown Logical; if TRUE, prints a markdown summary and the
#'   selected diagnostic rows. Returned values are not formatted.
#'
#' @return list(summary_tbl, projection_tbl)
#' @importFrom dplyr mutate select left_join group_by summarise add_row filter bind_rows arrange rename any_of
#' @importFrom tibble tibble
#' @importFrom purrr pmap_dfr map_dfr
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export
run_cancer_projection <- function(counts_5y,
                                  female_share, male_share,
                                  female_inc_wide, female_mort_5y,
                                  male_inc_wide,   male_mort_5y,
                                  study_start, study_end,
                                  age_min, age_max,
                                  begin_year = NULL, end_year,
                                  diag_years = NULL,
                                  print_markdown = TRUE) {

  # --- basic validations -----------------------------------------------------
  req_cols <- c("age_band", "N")
  if (!all(req_cols %in% names(counts_5y))) {
    stop("`counts_5y` must have columns: ", paste(req_cols, collapse = ", "))
  }
  rate_cols <- c("band", "rate_per100k")
  for (nm in c("female_inc_wide","female_mort_5y","male_inc_wide","male_mort_5y")) {
    x <- get(nm, inherits = FALSE)
    if (!all(rate_cols %in% names(x))) {
      stop("`", nm, "` must have columns: ", paste(rate_cols, collapse = ", "))
    }
  }
  if (!is.numeric(study_start) || !is.numeric(study_end) ||
      !is.numeric(age_min) || !is.numeric(age_max) ||
      !is.numeric(end_year)) {
    stop("Years and ages must be numeric.")
  }
  if (age_min >= age_max) stop("`age_min` must be < `age_max`.")
  if (is.null(begin_year)) begin_year <- floor((study_start + study_end) / 2)
  if (end_year < begin_year) stop("`end_year` must be >= `begin_year`.")

  # normalize sex shares
  total_share <- female_share + male_share
  if (total_share <= 0) stop("`female_share + male_share` must be > 0.")
  female_share <- female_share / total_share
  male_share   <- male_share   / total_share

  # diagnostics years
  if (is.null(diag_years)) {
    diag_years <- unique(c(begin_year, 2030, 2040))
    diag_years <- diag_years[diag_years >= begin_year & diag_years <= end_year]
  }

  ages_full <- seq.int(age_min, age_max)

  # --- strict validation of bands & coverage (catches 'o' vs '0') -----------
  # Cohort band strings: validate format only (coverage not required here)
  invisible(expand_bands_df(counts_5y, "age_band", ages_full))

  # Rate tables: validate strings, overlaps, and exact single-age coverage
  assert_bands_and_coverage(female_inc_wide, "band", ages_full, "Female incidence")
  assert_bands_and_coverage(female_mort_5y,  "band", ages_full, "Female mortality")
  assert_bands_and_coverage(male_inc_wide,   "band", ages_full, "Male incidence")
  assert_bands_and_coverage(male_mort_5y,    "band", ages_full, "Male mortality")

  # --- compute per sex ------------------------------------------------------
  female_out <- .project_expected_cases_one_sex(
    counts_5y, female_share, female_inc_wide, female_mort_5y,
    start_year = begin_year, end_year = end_year, ages_full = ages_full
  )
  male_out <- .project_expected_cases_one_sex(
    counts_5y, male_share,   male_inc_wide,   male_mort_5y,
    start_year = begin_year, end_year = end_year, ages_full = ages_full
  )

  # Summary (keep numeric; round only for display)
  tp <- paste0(begin_year, "-", end_year)
  summary_tbl <- tibble::tibble(
    Sex = c("Female", "Male", "Total"),
    `Time Period` = tp,
    `Expected Cases` = c(
      female_out$total_cases,
      male_out$total_cases,
      female_out$total_cases + male_out$total_cases
    )
  )

  # Diagnostics (selected years)
  diag_view <- dplyr::bind_rows(
    dplyr::mutate(female_out$diag_by_year, Sex = "Female"),
    dplyr::mutate(male_out$diag_by_year,   Sex = "Male")
  ) |>
    dplyr::select(
      Sex, year, alive_start, new_cases_year,
      deaths_year, aged_out_year, alive_end, cum_cases
    )

  diag_filtered <- diag_view |>
    dplyr::filter(.data$year %in% diag_years)

  if (isTRUE(print_markdown)) {
    .print_md_summary(summary_tbl)  # pretty print; keeps returned data numeric
    # Optional: integer-only display for diagnostics
    diag_display <- diag_filtered |>
      dplyr::mutate(
        alive_start    = round(alive_start),
        new_cases_year = round(new_cases_year),
        deaths_year    = round(deaths_year),
        aged_out_year  = round(aged_out_year),
        alive_end      = alive_start - deaths_year - aged_out_year,
        cum_cases      = cumsum(new_cases_year)
      )
    print(diag_display)
  }

  invisible(list(
    summary_tbl    = summary_tbl,
    projection_tbl = diag_filtered
  ))
}

# ==================== helpers: parsing/validation ============================

normalize_band <- function(b) {
  b <- gsub("\\s+", "", as.character(b))
  # normalize unicode dashes to ASCII hyphen
  b <- gsub("\u2013", "-", b, fixed = TRUE)  # en dash
  b <- gsub("\u2014", "-", b, fixed = TRUE)  # em dash
  b
}

valid_band_string <- function(b) {
  # allow only: NN-NN or NN+ (digits, hyphen, plus), no letters
  grepl("^[0-9]+(\\+|-[0-9]+)$", b)
}

parse_band_strict <- function(b, ages_full) {
  b0 <- normalize_band(b)
  if (!valid_band_string(b0)) {
    stop("Invalid age band '", b, "'. Use NN-NN or NN+ (digits only).")
  }
  if (grepl("\\+$", b0)) {
    lo <- as.numeric(sub("\\+.*$", "", b0))
    hi <- max(ages_full)
  } else {
    lo <- as.numeric(sub("-.*$", "", b0))
    hi <- as.numeric(sub("^.*-", "", b0))
  }
  if (is.na(lo) || is.na(hi) || hi < lo) stop("Failed to parse age band: '", b, "'.")
  c(lo = lo, hi = hi)
}

expand_bands_df <- function(tbl, col, ages_full) {
  raw <- tbl[[col]]
  if (is.null(raw)) stop("Column '", col, "' not found.")
  nb  <- normalize_band(raw)
  if (any(!valid_band_string(nb))) {
    bad <- unique(raw[!valid_band_string(nb)])
    stop("Found invalid age bands in '", col, "': ", paste(bad, collapse = ", "),
         ". Allowed: NN-NN or NN+; no letters (e.g., '0' not 'o').")
  }
  parsed <- t(vapply(nb, parse_band_strict, numeric(2), ages_full = ages_full))
  out <- data.frame(band_raw = raw, lo = parsed[, "lo"], hi = parsed[, "hi"], stringsAsFactors = FALSE)
  out[order(out$lo, out$hi), , drop = FALSE]
}

detect_overlaps <- function(bands_df) {
  # assumes sorted by lo, hi
  df <- bands_df[order(bands_df$lo, bands_df$hi), ]
  overlaps <- which(df$lo[-1] <= df$hi[-nrow(df)])
  overlaps
}

check_open_ended <- function(bands_df) {
  open_idx <- grep("\\+$", normalize_band(bands_df$band_raw))
  if (length(open_idx) > 1) {
    stop("Multiple open-ended bands found (e.g., 85+ and 90+). Keep exactly one.")
  }
  invisible(TRUE)
}

assert_exact_single_age_coverage <- function(bands_df, ages_full, label) {
  # Every age in ages_full must be covered by exactly one band
  cover_counts <- integer(length(ages_full))
  names(cover_counts) <- ages_full
  for (i in seq_len(nrow(bands_df))) {
    lo <- bands_df$lo[i]; hi <- bands_df$hi[i]
    idx <- which(ages_full >= lo & ages_full <= hi)
    if (length(idx)) cover_counts[idx] <- cover_counts[idx] + 1L
  }
  zero_ages  <- ages_full[cover_counts == 0L]
  multi_ages <- ages_full[cover_counts > 1L]

  if (length(zero_ages)) {
    stop(label, ": rate coverage is missing for ages ",
         paste(range(zero_ages), collapse = "-"),
         " (example ages: ", paste(head(zero_ages, 10), collapse = ", "), ").")
  }
  if (length(multi_ages)) {
    stop(label, ": overlapping rate bands detected around ages ",
         paste(range(multi_ages), collapse = "-"),
         ". Remove overlaps so each age maps to exactly one band.")
  }
  invisible(TRUE)
}

assert_bands_and_coverage <- function(rate_tbl, band_col, ages_full, label) {
  bands_df <- expand_bands_df(rate_tbl, band_col, ages_full)
  # explicit overlap check (fast fail, clearer message)
  ol <- detect_overlaps(bands_df)
  if (length(ol)) {
    bad_pairs <- paste0("['", bands_df$band_raw[ol], "' overlaps '", bands_df$band_raw[ol + 1], "']")
    stop(label, ": overlapping bands: ", paste(bad_pairs, collapse = ", "), ".")
  }
  check_open_ended(bands_df)
  assert_exact_single_age_coverage(bands_df, ages_full, label)
  invisible(TRUE)
}

# ==================== pretty printing (display rounding) =====================

.print_md_summary <- function(df) {
  stopifnot(is.data.frame(df))
  fmt <- df
  fmt$Sex <- ifelse(fmt$Sex == "Total", "**Total**", fmt$Sex)
  disp_vals <- round(df$`Expected Cases`)  # whole numbers for display
  fmt$`Expected Cases` <- ifelse(
    fmt$Sex == "**Total**",
    paste0("**", format(disp_vals, big.mark = ","), "**"),
    format(disp_vals, big.mark = ",")
  )

  cat("\n| Sex | Time Period | Expected Cases |\n")
  cat("|:---|:-----------|---------------:|\n")
  apply(fmt, 1, function(row) {
    cat("| ", paste(row[c("Sex","Time Period","Expected Cases")], collapse = " | "), " |\n", sep = "")
    invisible(NULL)
  })
}

# ==================== core expansion & projection ===========================

.expand_counts_to_single_years <- function(counts_5y, sex_prop, ages_full) {
  # Expand each 5-year band to single-year counts (equal split within the band)
  rows <- purrr::pmap_dfr(counts_5y, function(age_band, N) {
    rng <- parse_band_strict(age_band, ages_full)
    lo <- rng["lo"]; hi <- rng["hi"]
    w  <- hi - lo + 1
    tibble::tibble(age = lo:hi, n = (N * sex_prop) / w)
  })

  # Align to modeled range; fill missing ages with 0
  out <- tibble::tibble(age = ages_full) |>
    dplyr::left_join(
      rows |>
        dplyr::group_by(.data$age) |>
        dplyr::summarise(n = sum(.data$n), .groups = "drop"),
      by = "age"
    ) |>
    dplyr::mutate(alive_start = dplyr::coalesce(.data$n, 0)) |>
    dplyr::select(age, alive_start)

  # Conservation check
  total_expected <- sum(counts_5y$N) * sex_prop
  if (abs(sum(out$alive_start) - total_expected) > 1e-6) {
    stop("Expanded counts do not sum to the expected total after band expansion.")
  }

  out
}

.map_flat_rate_to_single_years <- function(rate_tbl, col_name, ages_full) {
  rt <- rate_tbl |>
    dplyr::rename(band = dplyr::any_of(c("age_band","band"))) |>
    dplyr::rename(rate_per100k = dplyr::any_of(c("inc_per100k","mort_per100k","rate_per100k")))

  parsed <- t(vapply(rt$band, parse_band_strict, numeric(2), ages_full = ages_full))
  bins <- rt |>
    dplyr::mutate(lo = parsed[, "lo"], hi = parsed[, "hi"]) |>
    dplyr::arrange(.data$lo, .data$hi)

  purrr::map_dfr(ages_full, function(a) {
    # vectorized match: which bands cover age 'a'?
    idx <- which(a >= bins$lo & a <= bins$hi)
    rate <- if (length(idx) > 0) bins$rate_per100k[idx[1]] else 0
    tibble::tibble(!!"age" := a, !!col_name := rate)
  })
}


.project_expected_cases_one_sex <- function(counts_5y, sex_prop, incidence_bins, mortality_bins,
                                            start_year, end_year, ages_full) {

  alive0   <- .expand_counts_to_single_years(counts_5y, sex_prop, ages_full)
  inc_age  <- .map_flat_rate_to_single_years(incidence_bins, "inc_per100k",  ages_full) |>
    dplyr::mutate(inc  = .data$inc_per100k  / 1e5)
  mort_age <- .map_flat_rate_to_single_years(mortality_bins, "mort_per100k", ages_full) |>
    dplyr::mutate(mort = .data$mort_per100k / 1e5)

  years <- seq.int(start_year, end_year)
  ages  <- ages_full
  age_min <- min(ages); age_max <- max(ages)

  alive    <- matrix(0, nrow = length(ages), ncol = length(years), dimnames = list(ages, years))
  exp_case <- matrix(0, nrow = length(ages), ncol = length(years), dimnames = list(ages, years))

  # baseline population at begin_year
  alive[, as.character(start_year)] <- alive0$alive_start[match(ages, alive0$age)]

  diag <- tibble::tibble(
    year = integer(), alive_start = double(), new_cases_year = double(),
    deaths_year = double(), aged_out_year = double(), alive_end = double(),
    cum_cases = double()
  )
  cum_cases <- 0

  for (t in seq_along(years)) {
    yr <- years[t]; ya <- as.character(yr)
    ia <- inc_age$inc[match(ages, inc_age$age)]
    ma <- mort_age$mort[match(ages, mort_age$age)]

    alive_start <- sum(alive[, ya])

    # expected new cases
    exp_case[, ya] <- alive[, ya] * ia
    new_cases_year <- sum(exp_case[, ya])

    # deaths
    deaths_vec  <- alive[, ya] * pmax(0, ma)
    deaths_year <- sum(deaths_vec)

    # survivors after deaths
    survivors <- pmax(alive[, ya] - deaths_vec, 0)

    # survivors at the top modeled age leave the range next year
    aged_out_year <- survivors[as.character(age_max)]

    # shift survivors to next year (aging +1), dropping the oldest
    if (t < length(years)) {
      next_yr <- years[t + 1]
      alive[as.character((age_min + 1):age_max), as.character(next_yr)] <-
        survivors[as.character(age_min:(age_max - 1))]
    }

    alive_end <- sum(survivors) - aged_out_year
    cum_cases <- cum_cases + new_cases_year

    diag <- dplyr::add_row(
      diag,
      year = yr,
      alive_start = alive_start,
      new_cases_year = new_cases_year,
      deaths_year = deaths_year,
      aged_out_year = aged_out_year,
      alive_end = alive_end,
      cum_cases = cum_cases
    )
  }

  list(
    total_cases = sum(exp_case),
    diag_by_year = diag
  )
}
