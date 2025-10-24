# Quiet R CMD check about NSE/tidyselect column names used in dplyr verbs
utils::globalVariables(c(
  "age", "alive_start", "Sex", "year",
  "new_cases_year", "deaths_year", "aged_out_year",
  "alive_end", "cum_cases"
))
