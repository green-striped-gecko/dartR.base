
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `dartR.base` <a href="https://green-striped-gecko.github.io/dartR/"><img src="man/figures/dartR_logo.png" align="right" height="140"/></a>

## An accessible genetic analysis platform for conservation, ecology and agriculture - dartR.base

<!-- badges: start -->
### Downloads
| package | downloads  |
|------------|------------|
|  dartR     | [![Downloads from Rstudio mirror](https://cranlogs.r-pkg.org/badges/grand-total/dartR)](https://cran.r-project.org/package=dartR) |
|  dartRverse  | [![Downloads from Rstudio mirror](https://cranlogs.r-pkg.org/badges/grand-total/dartRverse)](https://cran.r-project.org/package=dartRverse) |
| dartR.base   | [![Downloads from Rstudio mirror](https://cranlogs.r-pkg.org/badges/grand-total/dartR.base)](https://cran.r-project.org/package=dartR.base) |

### Repositories

| repo | status                                                                                                                                                                                                          |
|------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| main | [![main](https://github.com/green-striped-gecko/dartR.base/actions/workflows/check-standard.yml/badge.svg?branch=main)](https://github.com/green-striped-gecko/dartRverse/actions/workflows/check-standard.yml) |
| dev  | [![dev](https://github.com/green-striped-gecko/dartR.base/actions/workflows/check-standard.yml/badge.svg?branch=dev)](https://github.com/green-striped-gecko/dartR.base/actions/workflows/check-standard.yml)   |

Publication:
[![](https://img.shields.io/badge/doi-10.1111/1755--0998.12745-00cccc.svg)](https://doi.org/10.1111/1755-0998.12745)

Zenodo:
[![DOI](https://zenodo.org/badge/86640709.svg)](https://zenodo.org/badge/latestdoi/86640709)

<!-- badges: end -->

## Overview

`dartR.base` aims to support analysis of SNP data in a spatial context.
For example running a least-cost-path analysis or simple isolation by
distance plots. There is also a function that calculates spatial
autocorrelation.

dartR.base is one of the core packages of the dartRverse. Currently the
dartRverse consists of the following packages:

- dartR.base (report, filter and input/output functions, basic popgen)
- dartR.data
- dartR.sim (functions to simulate SNP data)
- dartR.spatial (spatial analysis)
- dartR.popgen (popgen analysis)
- dartR.sexlinked (identify and filter sexlinked markers)

`dartR` and its packages are a collaboration between the University of
Canberra, CSIRO, Diversity Arrays Technology, Arthur Rylah Institute for
Environmental Research and Monash University. `dartR` is supported with
funding from the ACT Priority Investment Program, CSIRO and the
University of Canberra.

<p align="center">
<img src='man/figures/institutions.png' width="800"/>
</p>

## Installation

For a normal install from CRAN use:

``` r
install.packages("dartR.base")
```

For hints and how to install github versions, check the help pages of
the dartRverse package.

## Contribute

If you want to help shape the future of `dartR`, [this
tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/Tutorial_0_dartR_for_the_Developer_2.0_19-Feb-22.pdf)
is for you.

## Citation

Please acknowledge `dartR` if you use it in your study. Copy and paste
the following code to the R console to retrieve the citation
information:

``` r
citation("dartR.base")
```

Check out our
[articles](https://github.com/green-striped-gecko/dartR/wiki/dartR-team-publications)
and our
[awards](https://github.com/green-striped-gecko/dartR/wiki/dartR-awards).

Have fun working with `dartR`!

Cheers,

Bernd, Arthur, Luis, Carlo, Diana & Olly
