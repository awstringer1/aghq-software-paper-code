# AGHQ Software Paper Code
Code for examples in AGHQ [software paper](https://arxiv.org/abs/2101.04468).

See also the `aghq` [package repository](https://github.com/awstringer1/aghq), and the `CRAN` [package page](https://cran.r-project.org/web/packages/aghq/index.html).

# Instructions

The following are instructions for reproducing the results in the paper.

## Install the package

First, install the `aghq` package stable version:

```
install.packages('aghq')
```

You will also need some other packages. Most can be obtained easily from `CRAN`, as follows:

- TMB
- tmbstan
- parallel
- glmmTMB
- geostatsp
- PrevMap
- geoR
- trustOptim
- numDeriv

There are three other packages that are more difficult to install. In all cases the script
checks for them and if they are not installed, throws a warning and skips the example in which
they are used:

- `ipoptr`: used in Example 4.2, Galactic Mass Estimation. Requires a working installation of `IPOPT`. See [here](https://coin-or.github.io/Ipopt/INSTALL.html).

- `INLA`: used as one method to compare to in Example 5.1, Loaloa. Not on `CRAN`, can be installed
from [this page](https://www.r-inla.org/download-install). 

- `EpiILMCT`: a dataset from this package is used 

Note that the `aghq` package itself has only a few dependencies and they are all on `CRAN`. But
some examples compare to other methods, use data from other packages, or, in one case, use one function from another
package in an intermediate step. 

## Run all examples at once

To run the examples I recommend you look at their individual files below. If you wish to run them all at once, you can do the following.

To run all the examples in the paper:

1. Install the latest version of `knitr`: `R -e 'install.packages("knitr")'`

2. Navigate to where you put the code and run the command `R -e 'knitr::spin("00-reproduce-all-results.R")'`.

This creates the files `00-reproduce-all-results.md` and `00-reproduce-all-results.html` which contain all the results from the paper.

The top of `00-reproduce-all-results.R` contains a command for installing all (but one) of the necessary packages, wrapped in a `if (FALSE)` statement so it won't be run on sourcing, but you can run it interactively if you need to.

The script will compile all the necessary `TMB` templates, which are including as files in the `aghq` package. It will create folders for each example inside the directory returned by `tempdir()`, and store all the plots, tables, and data (like `MCMC` results) from the paper there.

All the datasets used in the paper are available in the installed packages (including `aghq`).

## Run each example on its own

I recommend you look at the examples one at a time, in a clean `R` session.

Each code file creates results inside the directory returned by `tempdir()`. You can change this at the top of the file if you like.

All the files require packages:
  - `aghq`
  - `TMB`
  - `tmbstan`
  - The `parallel` package is used by `tmbstan` to run multiple chains. I haven't tested what happens on a Windows machine since I don't own one.

Some examples require additional packages as indicated in their individual scripts.

To consider each of the six examples in the paper individually, you can use the followig six code files, as follows:

1. Example 4.1 (infectious disease model): `02-disease.R`
2. Example 4.2 (galactic mass estimation): `01-astro.R`
3. Example 5.1 (Loaloa, without zero-inflation): `03-loaloa.R`
4. Example 5.2 (Loaloa, with zero-inflation): `05-loaloazip.R`
5. Example 2 (basic use): `07-basic-use.R`
6. Example 6 (Poisson ZIP from `glmmTMB`): `08-poisson-zip.R`

# Further Notes

The spatial examples require that the following lines run successfully:
```
cameroonBorderLL = raster::getData("GADM", country=c('CMR'), level=2)
nigeriaBorderLL = raster::getData("GADM", country=c('NGA'), level=2)
```
This depends on an external service and I have observed that this service is sometimes down. If you get a `cannot open connection` type error, check this. If it doesn't work you could just not plot the borders; the model fitting does not depend on this external service.
