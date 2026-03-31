#' Plot Mao's tau, empirical confidence intervals, and a fitted asymptotic model
#'
#' The \code{tau.plot()} function takes a species-by-sample incidence dataset
#' and visualizes Mao’s tau (the sample-based species-accumulation curve)
#' together with its empirical 95% confidence intervals. It also overlays a
#' fitted asymptotic model—either the Michaelis–Menten function
#' (\code{method = "MM"}) or a first-order rational function
#' (\code{method = "rational"})—and draws the corresponding asymptotic richness
#' estimate (limit) returned by \code{tau.limit()}. If the asymptotic limit
#' cannot be estimated, the fitted curve and horizontal limit line are omitted.
#'
#' For the empirical curve, the function computes Mao’s tau using
#' incidence-based rarefaction via \code{specaccum(method = "exact")}, optionally
#' accounting for the species occurrence structure when \code{conditioned = TRUE}.
#' Pointwise 95% confidence intervals are shown using \code{accum$sd} returned
#' by \code{specaccum()}. These intervals reflect sampling variation in the
#' order of samples and do not represent bootstrap uncertainty of the model
#' limit.
#'
#' @param data A species-by-sample incidence matrix or data frame.
#' @param method Character string specifying the asymptotic model to fit:
#'   \code{"MM"} for the Michaelis–Menten model or \code{"rational"} for the
#'   first-order rational function.
#' @param col Color used for points, confidence band, fitted curve, and limit.
#' @param new Logical; if \code{TRUE}, a new plot is created; if \code{FALSE},
#'   elements are added to an existing plot.
#' @param conditioned Logical; passed to \code{specaccum()}. If \code{TRUE},
#'   the algorithm accounts for the species occurrence structure when computing
#'   \code{accum$sd}.
#' @param ... Additional graphical parameters passed to plot() or lines().
#'
#' @return This function produces a graphical display and returns no value.
#'
#' @seealso \code{\link{tau.limit}}, \code{\link{tau.boot}}, \code{\link{specaccum}}
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
#' ## Plot two accumulation curves on the same graph
#' tau.plot(group1, method="rational", col="blue", new=TRUE, cond=TRUE)
#' tau.plot(group2, method="rational", col="red", new=FALSE, cond=TRUE)
#' }
#'
#' @export
#'
tau.plot <- function(data, method="MM", col="black", new=TRUE, conditioned=FALSE, ...){
  args<-list(...)
  mod<-tau.limit(data, method)
  lmt<-mod$limit
  accum <- suppressWarnings(
    accum<-specaccum(t(decostand(data, "pa")), method = "exact", conditioned=conditioned)
  )
  y<-accum$richness
  x<-seq_along(y)
  ymax <- if (is.finite(lmt)) lmt * 1.05 else max(y) * 1.1
  if (new){
    plot(x, y, col=col, pch=1, ylim=c(0, ymax), xlab="Number of samples", ylab="Number of species", ...)
  }else{
    points(x, y, col=col, pch=1)
  }
  upper<-y+1.96*accum$sd
  lower<-pmax(0, y-1.96*accum$sd)
  polygon(c(x, rev(x)), c(upper, rev(lower)), col=adjustcolor(col, alpha.f=0.15), border=NA)
  lines(x, upper, col=col, lty=3)
  lines(x, lower, col=col, lty=3)
  if (is.finite(lmt)) {
    if (!is.null(args$xlim)) x<-seq_len(args$xlim[2])
    if (method=="MM"){
      lines(x, mod$coef['a']*x/(mod$coef['b']+x), col=col, lty=1, ...)
    }else if (method=="rational"){
      lines(x, (mod$coef['a']+mod$coef['b']*x)/(1+mod$coef['c']*x), col=col, lty=1, ...)
    }
    abline(h=lmt, lty=5, col=col)
  }
}
