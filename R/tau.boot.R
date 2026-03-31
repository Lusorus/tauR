#' Bootstrap estimates of asymptotic species richness from Mao's tau curves
#'
#' Performs bootstrap resampling of samples in a species-by-sample dataset
#' and, for each resample, computes Mao’s tau values using incidence-based
#' rarefaction. The function then fits either the Michaelis–Menten model
#' (`method = "MM"`) or a first‑order rational function (`method = "rational"`)
#' to the resulting accumulation curve and returns a vector of bootstrap
#' estimates of asymptotic richness.
#'
#' @param data A species-by-sample matrix or data frame. Columns are treated
#'   as samples and rows as species. Abundances are converted to incidence
#'   (presence/absence) prior to analysis.
#' @param method Character string specifying the model used to approximate
#'   Mao’s tau values. Either `"MM"` for the Michaelis–Menten model or
#'   `"rational"` for a first-order rational function. Defaults to `"MM"`.
#' @param nperm Integer. Number of bootstrap resamples (with replacement)
#'   to generate. Defaults to 1000.
#'
#' @details
#' For each bootstrap resample of the sample set, Mao’s tau values are computed
#' using \code{specaccum(..., method = "exact")}. The selected model is then
#' fitted via nonlinear least squares. The asymptotic limit of the Michaelis–
#' Menten model is parameter a. For the rational model, the asymptotic
#' limit is obtained analytically as \eqn{b/c}, where \eqn{b} and \eqn{c} are
#' the fitted coefficients.
#'
#' @return A numeric vector of length \code{nperm} containing bootstrap
#'   estimates of asymptotic richness. Failed model fits return \code{NA}.
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
#' data(varespec)
#' data(varechem)
#'
#' ## Split samples by median humus depth
#' cutoff <- quantile(varechem$Humdepth, 0.5)
#' group1<-t(varespec[which(varechem$Humdepth <= cutoff), ])
#' group2<-t(varespec[which(varechem$Humdepth > cutoff), ])
#'
#' ## Bootstrap estimates of asymptotic richness
#' res1<-tau.boot(group1, "rational", nperm=200)
#' res2<-tau.boot(group2, "rational", nperm=200)
#'
#' ## One-sided test for differences in asymptotic richness between two groups
#'
#' if (mean(res1, na.rm = TRUE) >= mean(res2, na.rm = TRUE)) {
#'   delta <- na.omit(res1 - res2)
#' }else{
#'   delta <- na.omit(res2 - res1)
#' }
#'
#' alpha<-0.05
#' (q<-quantile(delta, alpha))
#'
#' ## If the lower alpha-quantile is > 0, the one-sided CI does not include 0
#' if (q>0) {
#'   print('The asymptotic richness differs between groups')
#' }else{
#'   print('No significant difference in asymptotic richness')
#' }
#'
#' ##Bootstrap p-value for the difference
#' mean(delta <= 0)
#' }
#'
#' @export
#'
tau.boot <- function (data, method = "MM", nperm = 1000){
data<-t(decostand(as.matrix(data), "pa"))
n <- nrow(data)
res<-c(rep(NA_real_, nperm))
for(i in seq_len(nperm)){
  xx<-data[sample(seq_len(n), n, replace=TRUE),]
  accum <- suppressWarnings(specaccum(xx, method = "exact"))
  y<-accum$richness
  x <- seq_along(y)
  if (method=="MM"){
    tryCatch({
      mod<-nls(y~(a*x)/(b+x), start=list(a=max(y), b=median(x)))
      res[i]<-coef(mod)[['a']]
    },error=function(msg){
      res[i]<-NA_real_
    })
  }else if (method=="rational"){
    tryCatch({
      mod<-nls(y~(a+b*x)/(1+c*x), start=list(a=0, b=1, c=0.1))
      limit <- coef(mod)[["b"]] / coef(mod)[["c"]]
      if (!is.finite(limit)) limit <- NA_real_
      res[i] <- limit
    },error=function(msg){
      res[i]<-NA_real_
    })
  } else {
    stop("Unknown method: use 'MM' or 'rational'")
  }
}
res
}
