# README

## Installation

### Prerequisites

You need a C++ compiler that sufficiently supports a C++11 standard, i.e.
GCC >= 5, or, recommended, LLVM >= 7; cf.

```
g++ --version
```

You need an `R` installation on your computer, then start R and run:

```
$ R
> install.packages(c('Rcpp', 'testthat', 'devtools', 'roxygen2'))
```

which should install all tools needed to build C/C++ based R packages.

*Note*: it's recommended to install above packages with suggested dependencies
(`dependencies=TRUE` argument of `install.packages`), however, beware that this
make take very long, in particular, for the `devtools` package.

### Build and install

To build and install the wrappers run in this folder:

```
make install
```

#### Troubleshooting

If you stumble upon problems while installing the `devtools` or `roxygen2` packages; or while installing `CASMAP` (e.g. `'no such file..'` error), then there are two additional steps you can follow

##### Simplified and explicit installation

In the subdirectory `r_wrapper` of the *root folder*, execute
```
make clean
make install_cwd    # build & install from a current working directory
```

##### Manual installation

In the subdirectory `r_wrapper` of the *root folder*, execute
```
make clean
```
to remove traces of a previous compilation

Build the package
```
Rscript -e "devtools::build()"
```
If the file `~/.R/Makevars` does not exist, create it. You may need to create the `~/.R` directory first.
Edit the file `~/.R/Makevars` and add the following contents:
```
CC=gcc-6
CXX=gcc-6
CXX1X=gcc-6
```
This allows R to compile the source code of the package using the `gcc-6` compiler. Refer to the environment variable `CC` that was set in the general installation steps.

Start R and install the package from source
```
$ R
> install.packages("../sigpatsearch_0.3.tar.gz", repos=NULL, type="source")
```

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
