# tauR

Tools for estimating total species richness, comparing expected species richness among habitats, and visualizing species accumulation curves.


## Installation

```r
# install.packages("remotes")
remotes::install_github("Lusorus/tauR")

# or install a specific archived version
remotes::install_url(
  "https://zenodo.org/record/19364998/files/Lusorus/tauR-v0.1.0.zip"
)
```


## Dependencies

The package **requires the `vegan` package** (Oksanen et al. 2012), which provides data structures and helper functions used in species accumulation analyses.  


## Quick start

```r
library(tauR)
?tauR

data(dune, package = "vegan")

## Asymptotic richness
tau.limit(t(dune))

## Basic CI manually
bt <- na.omit(tau.boot(t(dune)))
basic.CI <- c(
  2 * mean(bt) - quantile(bt, 0.025),
  mean(bt),
  2 * mean(bt) - quantile(bt, 0.975)
)
names(basic.CI) <- c("2.5%", "mean", "97.5%")
basic.CI

## Plot accumulation curve
tau.plot(t(dune))
```


## Functionality overview

- **tau.limit** — estimates the asymptotic limit of a species accumulation curve (Mao's tau) using either the Michaelis–Menten model or a first‑order rational function.

- **tau.boot** — generates bootstrap replicates of the asymptotic limit for confidence interval estimation and group comparisons.

- **tau.plot** — visualizes Mao’s tau, empirical confidence intervals, and fitted asymptotic models.

All functions expect community data matrices with species in rows and samples in columns.

## Citation

Seleznev, D. (2026). *Lusorus/tauR: tauR package*. Zenodo.  
[https://doi.org/10.5281/zenodo.19364997](https://doi.org/10.5281/zenodo.19364997)

This package is archived on Zenodo to ensure long‑term reproducibility. Please cite the **Concept DOI** (10.5281/zenodo.19364997) for general reference and the **Version DOI** (10.5281/zenodo.19364998) for exact reproducibility.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19364997.svg)](https://doi.org/10.5281/zenodo.19364997)


