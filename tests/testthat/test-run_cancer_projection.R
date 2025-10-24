# ------------------------------------------------------------
# test-run_cancer_projection.R
# ------------------------------------------------------------

# Load namespaces explicitly so tests are self-contained
library(testthat)
library(dplyr)
library(tibble)
library(expectedCases)

# --- Shared test fixturesult (tiny, fast, deterministic) -----------------------

counts_5y <- tibble::tibble(
  age_band = c("40-44","45-49","50-54","55-59"),
  N        = c(2694, 3480, 4166, 4964)
)

female_share <- 0.60
male_share   <- 0.40

female_inc_wide <- tibble::tibble(
  band = c("40-54","55-69","70-85"),
  rate_per100k = c(39, 102, 278)
)

female_mort_5y <- tibble::tibble(
  band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(87, 138, 214, 309, 455, 707, 1164, 2104, 4133, 13700)
)

male_inc_wide <- tibble::tibble(
  band = c("40-54","55-69","70-85"),
  rate_per100k = c(45, 148, 338)
)

male_mort_5y <- tibble::tibble(
  band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(155, 225, 337, 506, 789, 1190, 1888, 3239, 5974, 15415)
)

# --- Tests -------------------------------------------------------------------

test_that("run_cancer_projection returns expected structure", {
  result <- run_cancer_projection(
    counts_5y,
    female_share, male_share,
    female_inc_wide, female_mort_5y,
    male_inc_wide,   male_mort_5y,
    study_start = 2020, study_end = 2025,
    age_min = 40, age_max = 100,
    end_year = 2040,
    print_markdown = FALSE
  )

  # Structure
  expect_true(is.list(result))
  expect_true(all(c("summary_tbl", "projection_tbl") %in% names(result)))
  expect_s3_class(result$summary_tbl, "tbl_df")
  expect_s3_class(result$projection_tbl, "tbl_df")

  # Required columns in projection_tbl
  needed <- c("Sex","year","alive_start","new_cases_year","deaths_year",
              "aged_out_year","alive_end","cum_cases")
  expect_true(all(needed %in% names(result$projection_tbl)))
})

test_that("accounting identity holds at the yearly total level", {
  result <- run_cancer_projection(
    counts_5y,
    female_share, male_share,
    female_inc_wide, female_mort_5y,
    male_inc_wide,   male_mort_5y,
    study_start = 2020, study_end = 2025,
    age_min = 40, age_max = 100,
    end_year = 2040,
    print_markdown = FALSE
  )

  pt <- result$projection_tbl
  # alive_end + aged_out_year == alive_start − deaths_year

  lhs <- pt$alive_end + pt$aged_out_year
  rhs <- pt$alive_start - pt$deaths_year
  expect_equal(unname(lhs), unname(rhs), tolerance = 1e-8)
})

test_that("per-sex totals in summary equal sum of yearly new cases", {
  # Use a small, fully-covered scenario and include every year in diag_years
  counts_5y_local <- tibble::tibble(
    age_band = c("40-44","45-49"),
    N        = c(100, 100)
  )

  female_inc_wide_local <- tibble::tibble(band = c("40-44","45-49"), rate_per100k = c(50, 60))
  female_mort_5y_local  <- tibble::tibble(band = c("40-44","45-49"), rate_per100k = c(200, 250))
  male_inc_wide_local   <- tibble::tibble(band = c("40-44","45-49"), rate_per100k = c(50, 60))
  male_mort_5y_local    <- tibble::tibble(band = c("40-44","45-49"), rate_per100k = c(200, 250))

  begin <- 2020; finish <- 2025

  result <- run_cancer_projection(
    counts_5y_local,
    female_share = 0.6, male_share = 0.4,
    female_inc_wide_local, female_mort_5y_local,
    male_inc_wide_local,   male_mort_5y_local,
    study_start = begin, study_end = begin + 1,
    age_min = 40, age_max = 60,
    end_year = finish,
    diag_years = begin:finish,
    print_markdown = FALSE
  )

  # Sum of yearly new cases by sex across the full projection
  per_sex_sum <- result$projection_tbl |>
    dplyr::group_by(Sex) |>
    dplyr::summarise(total_new = round(sum(new_cases_year), 1), .groups = "drop")

  sum_tbl <- result$summary_tbl |>
    dplyr::filter(Sex %in% c("Female","Male","Total")) |>
    dplyr::mutate(`Expected Cases` = as.numeric(`Expected Cases`))


  # Compare Female / Male
  chk <- dplyr::left_join(
    per_sex_sum,
    dplyr::filter(sum_tbl, Sex %in% c("Female","Male")) |>
      dplyr::select(Sex, expected = `Expected Cases`),
    by = "Sex"
  )
  expect_equal(chk$total_new, chk$expected, tolerance = 1e-6)

  # And Total equals Female + Male
  total_expected <- sum_tbl$`Expected Cases`[sum_tbl$Sex == "Total"]
  expect_equal(total_expected, sum(chk$expected), tolerance = 1e-6)
})

test_that("open-ended rate bands like '85+' are handled and cumulative cases are non-decreasing", {
  counts_5y_local <- tibble::tibble(
    age_band = c("80-84","85-89"),
    N        = c(50, 30)
  )

  female_inc_wide_local <- tibble::tibble(band = c("80-84","85+"), rate_per100k = c(100, 200))
  female_mort_5y_local  <- tibble::tibble(band = c("80-84","85+"), rate_per100k = c(1000, 2000))
  male_inc_wide_local   <- tibble::tibble(band = c("80-84","85+"), rate_per100k = c(100, 200))
  male_mort_5y_local    <- tibble::tibble(band = c("80-84","85+"), rate_per100k = c(1000, 2000))

  result <- run_cancer_projection(
    counts_5y_local,
    female_share = 0.5, male_share = 0.5,
    female_inc_wide_local, female_mort_5y_local,
    male_inc_wide_local,   male_mort_5y_local,
    study_start = 2020, study_end = 2020,
    age_min = 80, age_max = 100,
    end_year = 2022,
    diag_years = 2020:2022,
    print_markdown = FALSE
  )

  # No NA in required columns
  need_cols <- c("Sex","year","alive_start","new_cases_year","deaths_year",
                 "aged_out_year","alive_end","cum_cases")
  expect_true(all(need_cols %in% names(result$projection_tbl)))
  expect_false(anyNA(result$projection_tbl[, need_cols]))

  # Cumulative cases are non-decreasing by year within each sex
  for (sx in unique(result$projection_tbl$Sex)) {
    v <- result$projection_tbl$cum_cases[result$projection_tbl$Sex == sx]
    expect_true(all(diff(v) >= -1e-12))
  }
})

