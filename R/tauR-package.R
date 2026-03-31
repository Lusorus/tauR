#' tauR: Tools for Asymptotic Richness Estimation and Accumulation Curves
#'
#' The package provides functions for estimating asymptotic species richness
#' using simple nonlinear models (Michaelis–Menten and first‑order rational),
#' computing bootstrap distributions of the asymptotic limit, and visualizing
#' sample‑based species‑accumulation curves (Mao’s tau).
#'
#' These tools support ecological and biodiversity analyses where
#' accumulation curves and nonparametric resampling are used to assess
#' asymptotic richness and compare communities.
#'
#' @details
#' The main user‑level functions are:
#'
#' \itemize{
#'   \item \code{\link{tau.limit}} — estimate the asymptotic limit of an
#'         accumulation curve using either the Michaelis–Menten or
#'         first‑order rational model.
#'
#'   \item \code{\link{tau.boot}} — obtain bootstrap replicates of the
#'         asymptotic limit for inference and group comparisons.
#'
#'   \item \code{\link{tau.plot}} — plot Mao’s tau, empirical confidence
#'         intervals, and fitted asymptotic models.
#' }
#'
#' @seealso
#' \code{\link{tau.limit}}, \code{\link{tau.boot}}, \code{\link{tau.plot}}
#'
#' @name tauR-package
#' @aliases tauR-package
#' @aliases tauR
#' @docType package
"_PACKAGE"

