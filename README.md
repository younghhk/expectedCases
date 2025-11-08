<!-- badges: start -->
[![R-CMD-check](https://github.com/younghhk/expectedCases/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/younghhk/expectedCases/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/younghhk/expectedCases/graph/badge.svg)](https://app.codecov.io/gh/younghhk/expectedCases)
[![MIT license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
<!-- badges: end -->
  
## expectedCases

Project expected cancer cases for a cohort using age-specific incidence and mortality rates. The package expands 5-year age-band counts to single-year ages, applies sex-specific rates, and simulates annual aging, incidence, and all-cause mortality.


📘 See the full vignette: [Tutorial: expectedCases](https://github.com/younghhk/expectedCases/blob/master/Vignettes_expectedCases.pdf)

## Installation

```{r}
# stable release (recommended)
remotes::install_github("younghhk/expectedCases@v0.1.2")

# or development head
remotes::install_github("younghhk/expectedCases")
```

## Overview

This quick start shows how to:

1.	prepare a cohort in 5-year age bands,  
	
2.	supply sex-specific incidence and mortality rate tables,
  
3.	run the projection, and
	
4.	interpret / visualize results.

### Requirements
	
  -  Cohort counts must be in 5-year bands (e.g., "40-44", "45-49", …).	
  -  Rate bands must be digits only (NN-NN or NN+), and must cover every single age in age_min:age_max exactly once (no gaps/overlaps).
  -  If you model up to, say, age_max = 100, include an open-ended incidence and mortality band like "85+" so older ages are covered.

## 1) Create the baseline cohort
```{r}
library(expectedCases)
library(tibble)
library(ggplot2)
library(dplyr)


counts_5y <- tibble(
  age_band = c("40-44","45-49","50-54","55-59"),
  N        = c(2694, 3480, 4166, 4964)
)
2) Define incidence and mortality rates
Incidence may be provided in broader bands (e.g., "40-54", "55-69"); the function maps these to single-year ages internally (flat within band). Mortality is typically in 5-year bands.
female_share <- 0.6
male_share   <- 0.4

# per 100k
female_inc_wide <- tibble(
  band = c("40-54", "55-69", "70-84", "85+"),
  rate_per100k = c(39, 102, 278, 400)
)

female_mort_5y <- tibble(
  band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(87, 138, 214, 309, 455, 707, 1164, 2104, 4133, 13700)
)

male_inc_wide <- tibble(
  band = c("40-54","55-69","70-84","85+"),
  rate_per100k = c(45, 148, 338, 500)
)

male_mort_5y <- tibble(
  band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(155, 225, 337, 506, 789, 1190, 1888, 3239, 5974, 15415)
)
```

## Validation behavior
- Band strings are strictly parsed: only NN-NN or NN+. Typos like "4o-44" (letter “o”) error out.
- Coverage is enforced: each single age in age_min:age_max must be covered by exactly one rate band (no gaps/overlaps). If you see a coverage error, add an open-ended band (e.g., "85+") or lower age_max.

## 3) Run the projection
```{r}
res <- run_cancer_projection(
  counts_5y,
  female_share = female_share, male_share = male_share,
  female_inc_wide, female_mort_5y,
  male_inc_wide,   male_mort_5y,
  study_start = 2020, study_end = 2025,
  age_min = 40, age_max = 100,
  end_year = 2040,
  diag_years = 2022:2040,   # return annual rows for inspection/plots
  print_markdown = FALSE
)

projection_tbl <- res$projection_tbl  # year-by-year details
summary_tbl    <- res$summary_tbl     # totals by sex & period
```

## Outputs
	•	summary_tbl: expected cases by Sex and time period (begin_year–end_year).
	•	projection_tbl: one row per Sex × Year with:
	◦	alive_start, new_cases_year, deaths_year, aged_out_year,
	◦	alive_end, and cum_cases (cumulative cases to that year).

## 4) Plot cumulative expected cases
```{r}
ggplot(projection_tbl, aes(year, cum_cases, color = Sex)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Expected Cumulative Cancer Cases", x = "Year", y = "Cumulative Expected Cases") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid.minor = element_blank())
```


## Tips & troubleshooting
- If you get a coverage error, add an "85+" tail (or reduce age_max) and ensure bands don’t overlap.



## Citation

If you use this code or method, please cite as:

Vogtmann, E. & Hong, G. (2025). Estimating Expected Colorectal Cancer Cases Using National Rates and Cohort Data. GitHub Repository.



## 🔒 Repository Access

This repository is currently restricted to authorized collaborators.

Last updated: October 2025

