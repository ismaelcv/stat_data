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

# define a custom conf_interval function
conf_interval <- function(x, y = NULL, conf_level = 0.95) {
  if (is.numeric(x)) {
     z <- stats::t.test(x, y, conf.level = conf_level)
     out <- as.numeric(z$conf.int)
     if (length(z$estimate) == 2) est <- z$estimate[1] - z$estimate[2] else est <- z$estimate
     out <- dplyr::data_frame(low = out[1], mean = est, high = out[2])
     return(out)
  } else if (inherits(x, "lm") | inherits(x, "lmrob")  ) {
    out <- stats::confint(x, level = conf_level)
    out <- dplyr::data_frame(variable = rownames(out), low = out[,1], mean = coef(x), high = out[,2])
    return(out)
  } else {
    stop("You must supply a numeric variable or an lm/lmrob object!")
  }
}

# define an r-squared wrapper for linear models
r_squared <- function(x) {
  if (inherits(x, "lm") | inherits(x, "lmrob")  ) {
    return(summary(x)$r.squared)
  } else if (inherits(x, "cluster")) {
    return(1 - x$tot.withinss / x$totss)
  } else {
     stop("You must supply an lm, lmrob, or cluster object!")
  }
}

# add predictions and residuals to the model
add_model <- function(data, model) {
  if (inherits(model, "lm") | inherits(model, "lmrob")  ) {

    data <- mutate(data, pred = predict(model, newdata = data),
                         resid = resid(model, newdata = data))
    data

  } else if (inherits(model, "cluster")) {



  } else {

     stop("You must supply an lm, lmrob, or cluster object!")

  }

  return(data)
}

# construct clustering model
clust_model <- function(data, num_clusters = 2) {
  data <- select_if(data, is.numeric)
  index <- which(names(data) %in% c("pred", "resid", "cluster", "tsne1", "tsne2"))
  if (length(index)) data <- data[,-index]
  data <- na.omit(data)
  data <- scale(data)
  out <- kmeans(as.matrix(data), centers = num_clusters, iter.max = 200L, nstart = 25L)
  class(out) <- c("cluster", class(out))
  out
}

print.cluster <- function(x, ...) {
  cat(sprintf("\nA cluster object with %d groups.\n", length(out$size)))
  cat(paste("  cluster sizes:", paste(out$size, collapse = ", ")))
  cat("\n\n")
}

# lm models
lm_basic <- function(data, formula) {
  out <- lm(formula, data)
  class(out) <- c("lm_basic", class(out))
  out
}

lm_robust <- function(data, formula) {
  out <- robustbase::lmrob(formula, data, method = "SMDM")
  class(out) <- c("lm_robust", class(out))
  out
}

lm_select <- function(data, formula) {
  mf <- lm(formula, data, method = "model.frame")
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)

  out <- cv.glmnet(x, y)
  class(out) <- c("lm_select", class(out))
  out
}

# lm models
bm_basic <- function(data, formula) {
  out <- glm(formula, data, family = binomial)
  class(out) <- c("bm_basic", class(out))
  out
}

bm_robust <- function(data, formula) {
  out <- robustbase::glmrob(formula, data, family = binomial)
  class(out) <- c("bm_robust", class(out))
  out
}

bm_select <- function(data, formula) {
  mf <- lm(formula, data, method = "model.frame")
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)

  out <- cv.glmnet(x, y, family = "binomial")
  class(out) <- c("bm_select", class(out))
  out
}

ca <- read_csv("https://statsmaths.github.io/stat_data/chicago_meta.csv")

lm_basic(ca, owner_ratio ~ perc_20_units + median_income)
lm_robust(ca, owner_ratio ~ perc_20_units + median_income)
lm_select(ca, owner_ratio ~ perc_20_units + median_income)

bm_basic(ca, owner_ratio > 43 ~ perc_20_units + median_income)
z <- bm_robust(ca, owner_ratio > 43 ~ perc_20_units + median_income)
bm_select(ca, owner_ratio > 43 ~ perc_20_units + median_income)







