# Estimating Expected  Cancer  Cases

This repository contains the methodology and `R` code used to estimate **expected cancer cases** in a cohort study, based on age-specific incidence and mortality rates derived from national statistics.

---

## Overview

The model estimates the number of new cancer cases expected to occur within a cohort during a user-defined time window (e.g., 2025–2040), accounting for competing mortality.  

At each single-year age `a` and year `t`:

- The number of **expected cancer cases** is calculated as  
  `new_cases[a, t] = survivors[a, t] * λ_c[a]`,  
  where `λ_c[a]` is the age-specific cancer incidence rate (per person-year).

- The number of **survivors** moving to the next year is  
  `survivors[a+1, t+1] = survivors[a, t] * (1 - λ_d[a])`,  
  where `λ_d[a]` is the all-cause mortality rate (per person-year).

This iterative process continues from `start_year` to `end_year`.  
The total expected cancer cases are then summed across all ages and years:

`E_total = Σ_a Σ_t (survivors[a, t] * λ_c[a])`


---

## Implementation in R


```r
# install.packages("devtools")  # if you don't have it yet
devtools::install_github("younghhk/expectedCases")
# load the package
library(expectedCases)
```

##  Create the baseline cohort

Each age band must represent a **5-year range** (e.g., “40–44”, “45–49”, etc.).
The function assumes equally spaced 5-year intervals for age expansion.  


```{r cohort}
counts_5y <- tibble::tibble(
  age_band = c("40-44","45-49","50-54","55-59"),
  N        = c(2694, 3480, 4166, 4964)
)
counts_5y
```

---

## Define incidence and mortality rates

Rates are specified separately for females and males.  
- **Incidence rates** are typically provided in *broader grouped age bands* (e.g., “0–24”, “25–39”, “40–54”, “55–69”, “70–85+”).  
- **Mortality rates** are usually specified in *detailed 5-year bands* (e.g., “40–44”, “45–49”, “50–54”, etc.).


```r
library(tibble)


# Sex proportions
female_share <- 0.6
male_share   <- 0.4

# Incidence and mortality rates per 100k (simplified integers)
female_inc_wide <- tibble(
  band = c("40-54", "55-69", "70-85"),
  rate_per100k = c(39, 102, 278)
)


female_mort_5y <- tibble(
  band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(87, 138, 214, 309, 455, 707, 1164, 2104, 4133, 13700)
)


male_inc_wide <- tibble(
  band = c("40-54","55-69","70-85"),
  rate_per100k = c(45, 148, 338)
)


male_mort_5y <- tibble(
  band = c("40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+"),
  rate_per100k = c(155, 225, 337, 506, 789, 1190, 1888, 3239, 5974, 15415)
)

```



```r
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
```

## Function Arguments

| Argument | Type | Description |
|:--|:--|:--|
| `counts_5y` | `tibble` | Baseline cohort counts by **5-year age bands** (e.g., "40-44", "45-49"), with one column `age_band` and one column `N` (population size). |
| `female_share` | `numeric` | Proportion of the total cohort that is female (between 0 and 1). |
| `male_share` | `numeric` | Proportion of the total cohort that is male (between 0 and 1). Must satisfy `female_share + male_share = 1`. |
| `female_inc_wide` | `tibble` | Female **incidence rates** per 100,000 population, specified in broader age bands (e.g., "40-54", "55-69", "70-85"). The function interpolates these to 5-year intervals internally. |
| `female_mort_5y` | `tibble` | Female **all-cause mortality rates** per 100,000 population, given in 5-year bands (e.g., "40-44", "45-49", …, "85+"). |
| `male_inc_wide` | `tibble` | Male incidence rates per 100,000 population, formatted like `female_inc_wide`. |
| `male_mort_5y` | `tibble` | Male all-cause mortality rates per 100,000 population, formatted like `female_mort_5y`. |
| `study_start` | `numeric` | Starting calendar year of the cohort (e.g., `2020`). |
| `study_end` | `numeric` | Ending year of data collection or baseline projection period (e.g., `2025`). Used to determine midpoint for internal timing. |
| `age_min` | `numeric` | Minimum age of the modeled cohort (default: `40`). |
| `age_max` | `numeric` | Maximum modeled age (default: `100`). Individuals exceeding this are dropped from the cohort. |
| `end_year` | `numeric` | Final projection year (e.g., `2040`). Determines length of simulation. |
| `begin_year` | `numeric` | Optional. Custom first projection year (default = midpoint of `[study_start, study_end]`). |
| `diag_years` | `numeric` vector | Optional. Specific years at which projections (e.g., 2022, 2030, 2040) are returned in the output. |
| `print_markdown` | `logical` | If `TRUE`, prints a formatted markdown summary table of expected cases in the console. Set `FALSE` for silent operation. |


### Return Value

The function returns a list with two tibbles:

| Output | Description |
|:--|:--|
| `summary_tbl` | Aggregated expected cancer cases by sex and projection period. |
| `projection_tbl` | Detailed year-by-year projection showing `alive_start`, `new_cases_year`, `deaths_year`, `aged_out_year`, `alive_end`, and `cum_cases`. |


### Example Output

The summary tibble shows projected **expected cancer cases** for a hypothetical cohort from **2022 to 2040**, by sex.

```r
result$summary_tbl

# A tibble: 3 × 3
  Sex     `Time Period` `Expected Cases`
  <chr>   <chr>                   <dbl>
1 Female  2022–2040               371.0
2 Male    2022–2040               306.3
3 Total   2022–2040               677.4
```

The projection tibble shows selected years with the number alive at the start of each year (`alive_start`), expected new cases (`new_cases_year`), deaths (`deaths_year`), individuals aging out of the modeled range (`aged_out_year`), survivors at year end (`alive_end`), and cumulative cases (`cum_cases`).

```r
result$projection_tbl

# A tibble: 6 × 8
  Sex     year alive_start new_cases_year deaths_year aged_out_year alive_end cum_cases
  <chr>  <int>       <dbl>          <dbl>       <dbl>         <dbl>     <dbl>     <dbl>
1 Female  2022      14486.           12.1        52.8             0    14434.      12.1
2 Female  2030      13890.           19.4       111.              0    13779.     141. 
3 Female  2040      12044.           23.8       301.              0    11743.     371. 
4 Male    2022       9658.           11.0        58.7             0     9599.      11.0
5 Male    2030       9018.           16.3       113.              0     8904.     123. 
6 Male    2040       7366.           18.3       237.              0     7129.     306. 
```



## Citation

If you use this code or method, please cite as:

Vogtmann, E. & Hong, G. (2025). Estimating Expected Colorectal Cancer Cases Using National Rates and Cohort Data (ABC Study). GitHub Repository.



## 🔒 Repository Access

This repository is currently restricted to authorized collaborators.

Last updated: October 2025

