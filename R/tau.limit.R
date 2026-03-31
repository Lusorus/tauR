#' Estimate asymptotic species richness from Mao's tau accumulation curves
#'
#' Computes incidence-based Mao’s tau values for a species-by-sample dataset
#' and fits either the Michaelis–Menten model (`method = "MM"`) or a
#' first-order rational function (`method = "rational"`) to the resulting
#' accumulation curve. The function returns the estimated asymptotic richness
#' (limit) together with the fitted model coefficients.
#'
#' @param data A species-by-sample matrix or data frame. Columns are treated
#'   as samples and rows as species. Abundances are converted to incidence
#'   (presence/absence) prior to analysis.
#' @param method Character string specifying the model used to approximate
#'   Mao’s tau values. Either `"MM"` for the Michaelis–Menten model or
#'   `"rational"` for a first-order rational function. Defaults to `"MM"`.
#'
#' @details
#' Mao’s tau values are computed using \code{specaccum(..., method = "exact")}.
#' Input must be species × samples; the function internally transposes to samples
#' × species. For \code{method = "MM"}, the model \eqn{(a x)/(b + x)} is fitted
#' via nonlinear least squares. The asymptotic limit of the Michaelis–Menten model
#' is parameter a. For \code{method = "rational"}, the model
#' \eqn{(a + b x)/(1 + c x)} is fitted, and the asymptotic limit is obtained
#' symbolically as \eqn{b/c}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{limit}{Estimated asymptotic richness (numeric).}
#'     \item{coef}{Named numeric vector of fitted model coefficients.}
#'   }
#'   If model fitting fails, both elements are returned as \code{NA}.
#'
#' @references
#' Colwell, R.K., et al. (2012). Models and estimators linking individual-based
#' and sample-based rarefaction, extrapolation and comparison of assemblages.
#' \emph{Journal of Plant Ecology}.
#'
#' Ratkowsky, D.A. (1990). \emph{Handbook of Nonlinear Regression Models}.
#'
#' Schoener, T.W. (1976). The species–area relation within archipelagos.
#'
#' @examples
#' \donttest{
#'
#' data(dune)
#'
#' ## Asymptotic richness using two models
#' tau.limit(t(dune), method = "MM")
#' tau.limit(t(dune), method = "rational")
#'
#' ## Percentile and Basic CI manually
#'
#' bt<-na.omit(tau.boot(t(dune), "rational", nperm=200))
#' alpha <- 0.05
#'
#' (perc.CI<-quantile(bt, probs = c(alpha/2, 0.5, 1 - alpha/2), na.rm = TRUE))
#'
#' basic.CI <- c(2*mean(bt)-perc.CI[3], mean(bt), 2*mean(bt)-perc.CI[1])
#' names(basic.CI)<-c("2.5%", "mean", "97.5%")
#' basic.CI
#' }
#'
#' @export

tau.limit <- function (data, method = "MM"){
  data<-t(decostand(as.matrix(data), "pa"))
  accum <- suppressWarnings(specaccum(data, method = "exact"))
  y<-accum$richness
  x<-seq_along(y)
  mod <- NULL
  if (method=="MM"){
    out<-tryCatch({
      mod<-nls(y~(a*x)/(b+x), start=list(a=max(y), b=median(x)))
      list(limit = coef(mod)[["a"]], coef = coef(mod))
    },error=function(msg){
      list(limit=NA_real_, coef=c(a=NA_real_, b=NA_real_))
    })
  }else if (method=="rational"){
    out<-tryCatch({
      mod<-nls(y~(a+b*x)/(1+c*x), start=list(a=0, b=1, c=0.1))
      limit<-coef(mod)[["b"]]/coef(mod)[["c"]]
      if (!is.finite(limit)) limit <- NA_real_
      list(limit=limit, coef = coef(mod))
    },error=function(msg){
      list(limit=NA_real_, coef=c(a = NA_real_, b = NA_real_, c = NA_real_))
    })
  }
  out
}

