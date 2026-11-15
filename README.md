
<!-- README.md is generated from README.Rmd. Please edit that file -->

# basketBH

<!-- badges: start -->

<!-- badges: end -->

## Overview

`basketBH` is an R package for the design and evaluation of basket
clinical trials based on Simon-type two-stage designs and the
Benjamini-Hochberg (BH) multiple-testing procedure.

The package implements exact calculations of key operating
characteristics, including:

- True positive rate (TPR)
- Expected number of true rejections
- Expected number of false rejections
- Expected sample size under the null
- Optimal and minimax basket-trial designs

The methodology is described in the accompanying manuscript:

> Azriel, Rinott and Posch. *Information Borrowing by the BH Procedure
> in Basket Trials*.

## Installation

### Install from GitHub

``` r
install.packages("remotes")

remotes::install_github("davidazriel/basketBH")
```

Replace `"davidazriel/basketBH"` with the actual repository name if it
differs.

### Install a local development version

From the package directory:

``` r
devtools::install()
```

## Quick Start

The main function in the package is `design_bh()`, which searches for
Simon-type two-stage basket-trial designs satisfying user-specified
power and false-discovery-rate requirements.

``` r
library(basketBH)

design_bh(
  K = 10,
  p0 = 0.15,
  pa = 0.35,
  target_alpha = 0.20,
  target_power = 0.70,
  min_n = 10,
  max_n = 40
)
```

The function returns:

- An optimal design (minimum expected sample size under the null)
- A minimax design (minimum maximum sample size)

## Example: Computing Basket-Level Power

``` r
bh_tpr(
  n1 = 11,
  r1 = 2,
  n2 = 19,
  p0 = 0.15,
  pa = 0.35,
  K = 10,
  alpha = 0.20,
  K1 = 1
)
```

This computes the basket-level true positive rate when exactly one
basket is active.

## Package Structure

The package currently provides functions for:

### Two-stage design calculations

``` r
pmf_two_stage()
expected_sample_size()
```

### Exact BH operating characteristics

``` r
exact_bh_expectations()
bh_tpr()
```

### Design search

``` r
design_bh()
```

## Reproducibility

The package is intended to support the numerical results presented in
the accompanying manuscript. Future releases will include vignettes
reproducing examples and simulation studies from the paper.

## Development Status

The package is currently under active development.

Planned additions include:

- Additional operating characteristics
- Simulation tools
- Graphical summaries
- Reproduction of manuscript examples and tables
- Additional basket-trial methodologies

## License

MIT License.
