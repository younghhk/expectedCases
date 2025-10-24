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
#'   for **female** incidence and mortality.
#' @param male_inc_wide,male_mort_5y Same for **male**.
#' @param study_start,study_end Integers, recruitment window.
#' @param age_min,age_max Integers defining modeled age range.
#' @param begin_year Optional integer; default = floor(midpoint of
#'   `study_start` and `study_end`).
#' @param end_year Final projection year (inclusive).
#' @param diag_years Optional integer vector of years to return in diagnostics;
#'   default includes `begin_year`, 2030, 2040 (kept within range).
#' @param print_markdown Logical; if `TRUE`, prints a markdown summary and the
#'   selected diagnostic rows. Returned values are not formatted.
#'
#' @return A list with:
#' \describe{
#'   \item{summary_tbl}{Tibble with columns `Sex`, `Time Period`,
#'         and **numeric** `Expected Cases`.}
#'   \item{diag_filtered}{Tibble of selected years with columns:
#'         `Sex`, `year`, `alive_start`, `new_cases_year`, `deaths_year`,
#'         `aged_out_year`, `alive_end`, `cum_cases`.}
#' }
#'
#' @examples
#' \dontrun{
#' # See README for a full runnable example
#' }
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

  # --- validations -----------------------------------------------------------
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

  # --- compute per sex ------------------------------------------------------
  female_out <- .project_expected_cases_one_sex(
    counts_5y, female_share, female_inc_wide, female_mort_5y,
    start_year = begin_year, end_year = end_year, ages_full = ages_full
  )
  male_out <- .project_expected_cases_one_sex(
    counts_5y, male_share,   male_inc_wide,   male_mort_5y,
    start_year = begin_year, end_year = end_year, ages_full = ages_full
  )

  # Numeric summary (keep numeric; format only when printing)
  tp <- paste0(begin_year, "\u2013", end_year)  # en dash
  summary_tbl <- tibble::tibble(
    Sex = c("Female", "Male", "Total"),
    `Time Period` = tp,
    `Expected Cases` = c(
      round(female_out$total_cases, 1),
      round(male_out$total_cases, 1),
      round(female_out$total_cases + male_out$total_cases, 1)
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
    print(diag_filtered)
  }

  invisible(list(
    summary_tbl   = summary_tbl,
    projection_tbl = diag_filtered
  ))
}

# ---- internal helpers -------------------------------------------------------

# Pretty markdown summary (formats only for display)
.print_md_summary <- function(df) {
  stopifnot(is.data.frame(df))
  fmt <- df
  fmt$Sex <- ifelse(fmt$Sex == "Total", "**Total**", fmt$Sex)
  fmt$`Expected Cases` <- ifelse(
    fmt$Sex == "**Total**",
    paste0("**", format(round(df$`Expected Cases`, 1), nsmall = 1), "**"),
    format(round(df$`Expected Cases`, 1), nsmall = 1)
  )

  cat("\n| Sex | Time Period | Expected Cases |\n")
  cat("|:---|:-----------|---------------:|\n")
  apply(fmt, 1, function(row) {
    cat("| ", paste(row[c("Sex","Time Period","Expected Cases")], collapse = " | "), " |\n", sep = "")
    invisible(NULL)
  })
}

.expand_counts_to_single_years <- function(counts_5y, sex_prop, ages_full) {
  # Parse a single band like "40-44", "85+", "40-44", " 40-44 "
  parse_band_5y <- function(b) {
    b <- gsub("\\s+", "", as.character(b))
    b <- gsub("-", "-", b)         # normalize en-dash
    if (grepl("\\+$", b)) {
      lo <- as.numeric(sub("\\+.*$", "", b))
      hi <- max(ages_full)
    } else if (grepl("^\\d+-\\d+$", b)) {
      lo <- as.numeric(sub("-.*$", "", b))
      hi <- as.numeric(sub("^.*-", "", b))
    } else {
      # fallback: extract numbers; if only one, treat as single-year band
      nums <- as.numeric(unlist(regmatches(b, gregexpr("\\d+", b))))
      if (length(nums) >= 2) { lo <- nums[1]; hi <- nums[2] }
      else if (length(nums) == 1) { lo <- nums[1]; hi <- nums[1] }
      else { lo <- NA_real_; hi <- NA_real_ }
    }
    c(lo = lo, hi = hi)
  }

  # Expand each 5y band to single-year counts (split evenly)
  rows <- purrr::pmap_dfr(counts_5y, function(age_band, N) {
    rng <- parse_band_5y(age_band)
    lo <- rng["lo"]; hi <- rng["hi"]



    # Defensive guard: if parse failed, keep zero rows for that band (no warning)
    if (is.na(lo) || is.na(hi) || hi < lo) {
      return(tibble::tibble(age = integer(0), n = numeric(0)))
    }

    w <- hi - lo + 1
    tibble::tibble(age = lo:hi, n = (N * sex_prop) / w)
  })

  # Align to the modeled age range; fill missing ages with 0
  tibble::tibble(age = ages_full) |>
    dplyr::left_join(
      rows |>
        dplyr::group_by(.data$age) |>
        dplyr::summarise(n = sum(.data$n), .groups = "drop"),
      by = "age"
    ) |>
  dplyr::mutate(alive_start = dplyr::coalesce(.data$n, 0)) |>
    dplyr::select(age, alive_start)
}


.map_flat_rate_to_single_years <- function(rate_tbl, col_name, ages_full) {
  rt <- rate_tbl |>
    dplyr::rename(band = dplyr::any_of(c("age_band","band"))) |>
    dplyr::rename(rate_per100k = dplyr::any_of(c("inc_per100k","mort_per100k","rate_per100k")))

  parse_band <- function(b) {
    b <- gsub("\\s+", "", as.character(b))
    b <- gsub("\u2013", "-", b)  # normalize en-dash if present
    if (grepl("\\+$", b)) {
      lo <- as.numeric(sub("\\+.*$", "", b)); hi <- max(ages_full)
    } else if (grepl("^\\d+-\\d+$", b)) {
      lo <- as.numeric(sub("-.*$", "", b))
      hi <- as.numeric(sub("^.*-", "", b))
    } else {
      nums <- as.numeric(unlist(regmatches(b, gregexpr("\\d+", b))))
      if (length(nums) >= 2) { lo <- nums[1]; hi <- nums[2] }
      else if (length(nums) == 1) { lo <- nums[1]; hi <- nums[1] }
      else { lo <- NA_real_; hi <- NA_real_ }
    }
    c(lo = lo, hi = hi)
  }

  parsed <- t(vapply(rt$band, parse_band, numeric(2)))
  bins <- rt |>
    dplyr::mutate(lo = parsed[, "lo"], hi = parsed[, "hi"]) |>
    dplyr::arrange(.data$lo, .data$hi)

  purrr::map_dfr(ages_full, function(a) {
    row <- dplyr::filter(bins, !is.na(.data$lo), !is.na(.data$hi), a >= .data$lo, a <= .data$hi)
    rate <- if (nrow(row) > 0) row$rate_per100k[1] else 0
    # dynamic column name without := / tidy-eval
    tibble::as_tibble(setNames(list(a, rate), c("age", col_name)))
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

