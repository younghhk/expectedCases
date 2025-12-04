<!-- badges: start -->
[![R-CMD-check](https://github.com/younghhk/expectedCases/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/younghhk/expectedCases/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/younghhk/expectedCases/graph/badge.svg)](https://app.codecov.io/gh/younghhk/expectedCases)
[![MIT license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
<!-- badges: end -->
  

[![Back to Hub](https://img.shields.io/badge/⬅️%20Back%20to%20Hub-2962FF?style=for-the-badge)](https://github.com/younghhk/NCI)


## expectedCases

This code projects expected cancer cases for a prospective cohort study using age-specific incidence and mortality rates. The package expands age-band counts (e.g., 30-34, 35-39, 40-44, etc years old) to single-year ages, applies sex-specific incidence and mortality rates, and simulates annual aging, cancer incidence, and all-cause mortality.

This repository is currently restricted to authorized collaborators until the program is fully verified.  


📘 See the full vignette: [Tutorial: expectedCases](https://github.com/younghhk/expectedCases/blob/master/vignette.pdf)

## Installation

```{r}
# stable release (recommended)
remotes::install_github("younghhk/expectedCases@v0.1.2")

# or development head
remotes::install_github("younghhk/expectedCases")
```

## Overview

This quick start shows how to:

1.	provide cohort age distribution in 5-year age bands,

2.	input sex distribution proportions for the cohort,  
	
3.	supply sex-specific incidence and mortality rate tables,
  
4.	run the projection, and
	
5.	interpret / visualize results.

### Requirements

- Cohort age distribution may use bands of any width (1-year, 5-year, 10-year, or mixed). Wider age bands may be less precise.

- Incidence and mortality rate bands must follow valid numeric formatting (NN-NN or NN+) and must cover every single age in age_min:age_max exactly once (no gaps or overlaps).

- Rates must be expressed per 100,000 population.

- If modeling to higher ages (e.g., age_max = 100), include an open-ended band such as "80+" so all ages above the last closed band are covered.

## 1) Create the baseline cohort
```{r}
library(expectedCases)
library(tibble)
library(ggplot2)
library(dplyr)

# example based on cohort data from the Cancer Prevention Study-3 Gut and Oral Microbiome Substudy (GOMS)
# PMID: 40889876
# input cohort age distribution
counts_5y <- tibble::tibble(
  age_band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74"),
  N        = c(1095, 1103, 1463, 2004, 2622, 2827, 2433)
)

# input sex distribution proportions for the cohort
female_share <- 0.769
male_share   <- 0.231

# input female cancer incidence rates (cases/100,000 population)
# incidence data from Globocan 2022 version 1.1 (https://gco.iarc.fr/today/en/dataviz/)
# colorectal cancer
female_inc_wide <- tibble(
  band = c("40-54","55-69","70+"),
  rate_per100k = c(41.6, 89.1, 143.8)
)

# input female all-cause mortality rates (cases/100,000 population)
# mortality data from WHO mortality database (https://www.who.int/data/data-collection-tools/who-mortality-database)
# US age- and sex-specific mortality rates for 2022
female_mort_5y <- tibble(
 band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(207.4,268.6,411.0,615.0,921.9,1291.2,1941.0,3176.6,5358.5,12841.0)
)

# input male cancer incidence rates (cases/100,000 population)
male_inc_wide <- tibble(
  band = c("40-54","55-69","70+"),
  rate_per100k = c(51.4, 119.6, 174.3)
)

# input male all-cause mortality rates (cases/100,000 population)
male_mort_5y <- tibble(
  band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(382.5,465.5,684.5,990.6,1479.5,2045.3,2897.5,4532.0,7103.3,13876.6)
)
```

## 3) Run the projection
```{r}
res <- run_cancer_projection(
  counts_5y,
  female_share = female_share, male_share = male_share,
  female_inc_wide, female_mort_5y,
  male_inc_wide,   male_mort_5y,
  study_start = 2020, study_end = 2023, # input the years for the start and end of study recruitment
  age_min = 40, age_max = 100,
  end_year = 2040,
  diag_years = 2022:2040,   # return diagnostic rows for these specific years
  print_markdown = FALSE
)

projection_tbl <- res$projection_tbl  # year-by-year details
summary_tbl    <- res$summary_tbl     # totals by sex & period
```

## Outputs

- summary_tbl: expected cases by Sex and time period (begin_year–end_year).
- projection_tbl: one row per Sex × Year with:
	- alive_start: number of individuals alive at the start of the time period
 	- new_cases_year: number of new cancer cases for the year at the end of the time period
  	- deaths_year: number of deaths of the year at the end of the time period
  	- aged_out_year: number of individuals who age out of the cohort (i.e., exceed the maximum age)
	- alive_end: number of individuals alive at the end of the time period
 	- cum_cases: cumulative number of cancer cases at the end of the time period

## 4) Plot cumulative expected cases
```{r}
ggplot(projection_tbl, aes(year, cum_cases, color = Sex)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Expected Cumulative Cancer Cases", x = "Year", y = "Cumulative Expected Cases") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid.minor = element_blank())
```


## Citation

If you use this code or method, please cite as:

Vogtmann, E. & Hong, G. H. (2025). Projecting Expected  Cancer Cases for Cohort Data. GitHub Repository.



## 🔒 Repository Access

Access is currently limited to authorized collaborators pending full program verification.

Last updated: December 2025

