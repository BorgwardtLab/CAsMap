# README

## Installation

### Prerequisites

You need a C++ compiler that sufficiently supports a C++11 standard, i.e.
GCC >= 5, or, recommended, LLVM >= 7; cf.

    $ g++ --version

You need an `R` installation on your computer, then run:

    $ R
    > install.packages(c('Rcpp', 'testthat', 'devtools', 'roxygen2'))

which should install all tools needed to build C/C++ based R packages.

*Note*: it's recommended to install above packages with suggested dependencies
(`dependencies=TRUE` argument of `install.packages`), however, beware that this
make take very long, in particular, for the `devtools` package.

### Build, install and test

To build and install the wrappers run in this folder:

    $ make install

To run then tests:

    $ make test

**Beware**: tests check output which depends on using `rand`/`srand` from
`stdlib` (WY methods), so despite setting same seed you might still get
_slightly_ different results between different versions of `stdlib` (different
OSs in particular). Checks should always pass for non-WY methods and with
C libraries provided by Apple LLVM version 8.0.0.

#### Troubleshooting

If you stumble upon problems while installing `devtools`, `roxygen2` packages,
or while installing this package (e.g. `'no such file..'` error), then try
simplified and explicit installation and tests calls:

    $ make clean
    $ make install_cwd # build & install from a current working directory
    $ (cd tests/ && Rscript testthat.R) # run tests w/o devtools

### Build with CRAN flgs and check

To just build the wrapper with extra CRAN curator's flags, which print more
warnings, run in this folder:

    $ make wrapper

To just package the wrapper run in this folder:

    $ make package

This will create a zipped TAR file `../sigpatsearch_PKGVER.tar.gz`,
where `PKGVER` is the current package version. You can check this package for
CRAN publication by running:

    $ make check

This creates the `sigpatsearch.Rcheck/` folder which contains also a R
vignette PDF file.

### Cleanup

To remove build and packaging files run:

    $ make clean

To uninstall the package from your system run:

    $ make uninstall
