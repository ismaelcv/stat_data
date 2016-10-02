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
check_and_load("tsne")
check_and_load("GGally")

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
  set.seed(1)

  data <- select_if(data, is.numeric)

  if (missing(formula)) {
    formula <- ~.
  }

  if (num_clusters < 2) stop("Must select at least two clusters")
  num_clusters <- as.integer(num_clusters[[1]])

  mf <- stats::lm(formula, data, method = "model.frame")
  mt <- attr(mf, "terms")
  xlevels <- .getXlevels(mt, mf)
  mm <- model.matrix(formula, data)
  X <- scale(mm)
  center <- attr(X,"scaled:center")
  scale <- attr(X,"scaled:scale")
  X <- X[,colnames(X) != "(Intercept)"]


  model <- flexclust::stepFlexclust(X, k = num_clusters, nrep = 10L, simple = FALSE)
  model@call <- match.call()
  model_output <- list(terms = mt, contrasts = NULL,
                             xlevels = xlevels, center = center,
                             scale = scale)

  Xdf <- dplyr::as_data_frame(X)
  Xdf$CLUSTER <- model@cluster
  Xdf <- group_by_(Xdf, "CLUSTER")
  z <- reshape2::melt(summarize_all(Xdf, mean), id.vars = c("CLUSTER"))

  attributes(model)$model_output <- model_output
  attributes(model)$formula <- formula
  attributes(model)$clsummary <- z

  model
}

.cm_predict <- function(x, newdata, ...) {
  if (missing(newdata))
    stop("You must supply a new data object in order to run prediction")

  # construct new data matrix
  object <- attributes(x)$model_output
  tt <- terms(object)
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, na.action = na.pass, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, mf)
  X <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)

  # scale the data matrix
  X <- scale(X, center = object$center, scale = object$scale)
  X <- X[,colnames(X) != "(Intercept)"]

  # apply the prediction algorithm
  predictor <- predict(x, X)

  predictor
}


.cm_resid <- function(x, newdata, response = TRUE, ...) {
  if (missing(newdata))
    stop("You must supply a new data object in order to run prediction")

  # construct new data matrix
  object <- attributes(x)$model_output
  tt <- terms(object)
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, na.action = na.pass, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, mf)
  X <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)

  # scale the data matrix
  X <- scale(X, center = object$center, scale = object$scale)
  X <- X[,colnames(X) != "(Intercept)"]

  # apply the prediction algorithm
  predictor <- predict(x, X)

  # get distance to the cluster center
  dists <- cbind(distEuclidean(X, x@centers), predictor)
  resid <- apply(dists, 1, function(v) v[v[length(v)]] )

  # add variation of the response, if needed
  if (response) {
    raw_resid <- sqrt(apply(X^2, 1, sum))
    attributes(resid)$response <- var(raw_resid, na.rm=TRUE)
  }

  # return the residual
  resid
}

add_clusters <- function(data, model) {

  predictor <- .cm_predict(model, data)
  mutate(data, cluster = predictor)

}

add_tsne <- function(data, model) {

  formula <- attributes(model)$formula

  mf <- stats::lm(formula, data, method = "model.frame")
  mt <- attr(mf, "terms")
  x <- stats::model.matrix(mt, mf)

  out <- tsne::tsne(x)
  dplyr::mutate(data, tsne1 = out[,1], tsne2 = out[,2])

}

explained_variance <- function(data, object){

  resid  <- .cm_resid(object, data, response = TRUE)
  raw <- attributes(resid)$response
  ev <- 1 - var(resid, na.rm = TRUE) / raw
  return(ev)

}

plot_clusters <- function(model) {

  z <- attributes(model)$clsummary
  ggplot(z, aes(variable, value)) +
    geom_bar(stat = "identity", position = "identity", show.legend = FALSE) +
    coord_flip() +
    facet_wrap( ~ CLUSTER) +
    xlab("") +
    ylab("") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

}

lm_basic <- function(formula, data) {
  out <- stats::lm(formula, data)
  out$call <- match.call()
  class(out) = c("lm_basic", "lm")
  out
}

lm_robust <- function(formula, data) {
  out <- robustbase::lmrob(formula, data, method = "MM")
  out$call <- match.call()
  class(out) = c("lm_robust", "lm")
  out
}

add_prediction <- function(data, object){
  pred <- predict(object, newdata = data)
  resid <- resid(object, newdata = data)

  model_name <- deparse(substitute(object))

  data[,paste(model_name,"pred",sep="_")] <- pred
  data[,paste(model_name,"resid",sep="_")] <- resid

  data
}

summary.lm_basic <- function (object, correlation = FALSE, symbolic.cor = FALSE,
    level = 0.95, mu = 0, ...)
{
    z <- object

    # MODIFIED
    ci_vals <- confint(object, level = level)

    p <- z$rank
    rdf <- z$df.residual
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        }
        else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/rdf
        ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))
        ans$residuals <- r
        ans$df <- c(0L, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0L, 4L)
        dimnames(ans$coefficients) <- list(NULL, c("Estimate",
            "Std. Error", "t value", "Pr(>|t|)"))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        return(ans)
    }
    if (is.null(z$terms))
        stop("invalid 'lm' object:  no 'terms' component")
    if (!inherits(object, "lm"))
        warning("calling summary.lm(<fake-lm-object>) ...")
    Qr <- stats:::qr.lm(object)
    n <- NROW(Qr$qr)
    if (is.na(z$df.residual) || n - p != z$df.residual)
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept"))
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) *
        1e-30)
        warning("essentially perfect fit: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- (est - mu)/se
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$residuals <- r
    ans$coefficients <- cbind(est,ci_vals, 2 * pt(abs(tval),
        rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]],
        c("Estimate", colnames(ci_vals), "Pr(>|t|)"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept"))
            1L
        else 0L
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
            df.int)/rdf)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
            numdf = p - df.int, dendf = rdf)
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,
        1)]
    if (correlation) {
        ans$correlation <- (R * resvar)/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
        ans$symbolic.cor <- symbolic.cor
    }
    if (!is.null(z$na.action))
        ans$na.action <- z$na.action
    class(ans) <- "summary.lm"
    ans
}

summary.lm_robust <- function (object, correlation = FALSE, symbolic.cor = FALSE,
    level = 0.95, mu = 0, ...)
{
    # MODIFIED
    ci_vals <- confint(object, level = level)

    if (is.null(object$terms))
        stop("invalid 'lmrob' object:  no terms component")
    p <- object$rank
    df <- object$df.residual
    sigma <- object[["scale"]]
    aliased <- is.na(coef(object))
    cf.nms <- c("Estimate", colnames(ci_vals), "Pr(>|t|)")
    if (p > 0) {
        n <- p + df
        p1 <- seq_len(p)
        se <- sqrt(if (length(object$cov) == 1L) object$cov else diag(object$cov))
        est <- object$coefficients[object$qr$pivot[p1]]
        tval <- (est - mu)/se
        ans <- object[c("call", "terms", "residuals", "scale",
            "rweights", "converged", "iter", "control")]
        if (!is.null(ans$weights))
            ans$residuals <- ans$residuals * sqrt(object$weights)
        ans$df <- c(p, df, NCOL(object$qr$qr))
        ans$coefficients <- if (ans$converged)
            cbind(est, ci_vals, 2 * pt(abs(tval), df, lower.tail = FALSE))
        else cbind(est, if (sigma <= 0)
            0
        else NA, NA, NA)
        dimnames(ans$coefficients) <- list(names(est), cf.nms)
        if (p != attr(ans$terms, "intercept")) {
            df.int <- if (attr(ans$terms, "intercept"))
                1L
            else 0L
            resid <- object$residuals
            pred <- object$fitted.values
            resp <- if (is.null(object[["y"]]))
                pred + resid
            else object$y
            wgt <- object$rweights
            ctrl <- object$control
            c.psi <- ctrl$tuning.psi
            psi <- ctrl$psi
            correc <- if (psi == "ggw") {
                if (isTRUE(all.equal(c.psi, c(-0.5, 1, 0.95,
                  NA))))
                  1.121708
                else if (isTRUE(all.equal(c.psi, c(-0.5, 1.5,
                  0.95, NA))))
                  1.163192
                else if (isTRUE(all.equal(c.psi, c(-0.5, 1, 0.85,
                  NA))))
                  1.33517
                else if (isTRUE(all.equal(c.psi, c(-0.5, 1.5,
                  0.85, NA))))
                  1.395828
                else lmrob.E(wgt(r), ctrl)/lmrob.E(r * psi(r),
                  ctrl)
            }
            else if (any(psi == robustbase:::.Mpsi.R.names) && isTRUE(all.equal(c.psi,
                .Mpsi.tuning.default(psi)))) {
                switch(psi, bisquare = 1.207617, welsh = 1.224617,
                  optimal = 1.068939, hampel = 1.166891, lqq = 1.159232,
                  stop("unsupported psi function -- should not happen"))
            }
            else lmrob.E(wgt(r), ctrl)/lmrob.E(r * psi(r), ctrl)
            resp.mean <- if (df.int == 1L)
                sum(wgt * resp)/sum(wgt)
            else 0
            yMy <- sum(wgt * (resp - resp.mean)^2)
            rMr <- sum(wgt * resid^2)
            ans$r.squared <- r2correc <- (yMy - rMr)/(yMy + rMr *
                (correc - 1))
            ans$adj.r.squared <- 1 - (1 - r2correc) * ((n - df.int)/df)
        }
        else ans$r.squared <- ans$adj.r.squared <- 0
        ans$cov <- object$cov
        if (length(object$cov) > 1L)
            dimnames(ans$cov) <- dimnames(ans$coefficients)[c(1,
                1)]
        if (correlation) {
            ans$correlation <- ans$cov/outer(se, se)
            ans$symbolic.cor <- symbolic.cor
        }
    }
    else {
        ans <- object
        ans$df <- c(0L, df, length(aliased))
        ans$coefficients <- matrix(NA, 0L, 4L, dimnames = list(NULL,
            cf.nms))
        ans$r.squared <- ans$adj.r.squared <- 0
        ans$cov <- object$cov
    }
    ans$aliased <- aliased
    ans$sigma <- sigma
    if (is.function(ans$control$eps.outlier))
        ans$control$eps.outlier <- ans$control$eps.outlier(nobs(object))
    if (is.function(ans$control$eps.x))
        ans$control$eps.x <- if (!is.null(o.x <- object[["x"]]))
            ans$control$eps.x(max(abs(o.x)))
    structure(ans, class = "summary.lm_robust")
}

print.summary.lm_robust <-
function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), showAlgo = TRUE,
    ...)
{
    cat("\nCall:\n", paste(deparse(x$call, width.cutoff = 72),
        sep = "\n", collapse = "\n"), "\n", sep = "")
    control <- robustbase:::lmrob.control.neededOnly(x$control)
    cat(" \\--> method = \"", control$method, "\"\n", sep = "")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    cat(if (!is.null(x$weights) && diff(range(x$weights)))
        "Weighted ", "Residuals:\n", sep = "")
    if (rdf > 5L) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (NCOL(resid) > 1)
            structure(apply(t(resid), 1, quantile), dimnames = list(nam,
                dimnames(resid)[[2]]))
        else setNames(quantile(resid), nam)
        print(rq, digits = digits, ...)
    }
    else print(resid, digits = digits, ...)
    if (length(x$aliased)) {
        if (!(x$converged)) {
            if (x$scale == 0) {
                cat("\nExact fit detected\n\nCoefficients:\n")
            }
            else {
                cat("\nAlgorithm did not converge\n")
                if (control$method == "S")
                  cat("\nCoefficients of the *initial* S-estimator:\n")
                else cat(sprintf("\nCoefficients of the %s-estimator:\n",
                  control$method))
            }
            printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
                ...)
        }
        else {
            if (nsingular <- df[3L] - df[1L])
                cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
                  sep = "")
            else cat("\nCoefficients:\n")
            coefs <- x$coefficients
            if (!is.null(aliased <- x$aliased) && any(aliased)) {
                cn <- names(aliased)
                coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                  colnames(coefs)))
                coefs[!aliased, ] <- x$coefficients
            }
            printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                na.print = "NA", ...)
            cat("\nRobust residual standard error:", format(signif(x$scale,
                digits)), "\n")
            if (!is.null(x$r.squared) && x$df[1] != attr(x$terms,
                "intercept")) {
                cat("Multiple R-squared: ", formatC(x$r.squared,
                  digits = digits))
                cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared,
                  digits = digits), "\n")
            }
            correl <- x$correlation
            if (!is.null(correl)) {
                p <- NCOL(correl)
                if (p > 1) {
                  cat("\nCorrelation of Coefficients:\n")
                  if (is.logical(symbolic.cor) && symbolic.cor) {
                    print(symnum(correl), abbr.colnames = NULL)
                  }
                  else {
                    correl <- format(round(correl, 2), nsmall = 2,
                      digits = digits)
                    correl[!lower.tri(correl)] <- ""
                    print(correl[-1, -p, drop = FALSE], quote = FALSE)
                  }
                }
            }
            cat("Convergence in", x$iter, "IRWLS iterations\n")
        }
        cat("\n")
    }
    else cat("\nNo Coefficients\n")
    invisible(x)
}
