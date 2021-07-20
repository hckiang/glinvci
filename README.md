# glinvci
[![Build Status](https://travis-ci.com/hckiang/glinvci.svg?branch=main)](https://travis-ci.com/hckiang/glinvci)

`glinv` is an R package which enables computation of asymptotic confidence intervals for a class of continuous-time
Gaussian branching processes along known phylogenies, including the OU process commonly used in phylogenetic
comparative methods.

## Installation

```{r}
devtools::install_github("hckiang/glinvci")
```

## Usage

The example section of `?glinvci` contains the most basic usage in unrestricted OU process model.
Handling of missing data is documented in `?ou_haltlost`. To restrict parameters, take a look
at the examples in `?ou_diagH`.
