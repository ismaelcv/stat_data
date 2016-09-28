# install packages if not already present
check_and_load <- function(pkg) {
  z <- capture.output( if (!require(pkg, character.only = TRUE)) install.packages(pkg), type = "message" )
  z <- capture.output( library(pkg, character.only = TRUE), type = "message")
}

check_and_load("ggplot2")
check_and_load("dplyr")
check_and_load("readr")
check_and_load("ggrepel")
check_and_load("robustbase")
check_and_load("flexclust")
check_and_load("glmnet")

# define a custom group_summary function
group_summarize <- function(.data, ...) {

  data <- group_by_(.data, .dots = lazyeval::lazy_dots(...))
  group_vars <- sapply(attributes(data)$vars, as.character)

  data <- select_if(data, is.numeric)
  results <- summarize_all(data, funs(mean, median, sd, sum))

  if (ncol(data) <= length(group_vars) + 1) {
    these <- (names(results) %in% group_vars)
    names(results)[!these] <- paste0(setdiff(names(data),group_vars), "_", names(results)[!these])
  }

  results$n <- summarize(data, n = n())$n
  ungroup(results)
}

cm_basic <- function(data, formula, num_clusters = 2) {
  out <- stats::lm(formula, data, method = "model.frame")
  mm <- model.matrix(formula, data)
  X <- scale(mm)

  model <- flexclust::stepFlexclust(X, k = num_clusters, nrep = 10L, simple = TRUE)
  model$call <- match.call()
  model$terms <- out$terms
  model$contrasts <- out$contrasts
  model$xlevels <- out$xlevels
  model$center <- attr(X,"scaled:center")
  model$scale <- attr(X,"scaled:scale")
  class(model) <- c("cm_basic", "kccasimple")

  model
}




ca <- read_csv("https://statsmaths.github.io/stat_data/chicago_meta.csv")

lm_basic(ca, owner_ratio ~ perc_20_units + median_income)
lm_robust(ca, owner_ratio ~ perc_20_units + median_income)
lm_select(ca, owner_ratio ~ perc_20_units + median_income)

bm_basic(ca, owner_ratio > 43 ~ perc_20_units + median_income)
z <- bm_robust(ca, owner_ratio > 43 ~ perc_20_units + median_income)
bm_select(ca, owner_ratio > 43 ~ perc_20_units + median_income)







