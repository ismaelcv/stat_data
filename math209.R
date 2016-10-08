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
check_and_load("ggmap")

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

lm_robust <- function(formula, data, ...) {
  out <- robustbase::lmrob(formula, data, setting = "KS2014")
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


summary.lm_robust <- function (object, correlation = FALSE, symbolic.cor = FALSE,
    conf.level = 0.95, mu = 0, ...)
{
    # MODIFIED
    ci_vals <- confint(object, level = conf.level)

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

    # ADJUSTMENTS
    #ans$coefficients[,4] <- p.adjust(ans$coefficients[,4])

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
            printCoefmat(x$coef, digits = digits, signif.stars = FALSE,
                has.Pvalue = TRUE,
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
            printCoefmat(coefs, digits = digits, signif.stars = FALSE,
                na.print = "NA", has.Pvalue = TRUE, ...)
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


