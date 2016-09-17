# install packages if not already present
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("readr")) install.packages("readr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggrepel")) install.packages("ggrepel")

# load the packages
library(MASS)
library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)

# define a custom group_summary function
group_summarize <- function(.data, ...) {

  data <- group_by_(.data, .dots = lazyeval::lazy_dots(...))
  group_vars <- sapply(attributes(z)$vars, as.character)

  data <- select_if(data, is.numeric)
  results <- summarize_all(data, funs(mean, median, sd, sum))

  if (ncol(data) <= length(group_vars) + 1) {
    these <- (names(results) %in% group_vars)
    names(results)[!these] <- paste0(setdiff(names(data),group_vars), "_", names(results)[!these])
  }

  results$n <- summarize(data, n = n())$n
  ungroup(results)
}

