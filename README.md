# CASMAP
**C**ombinatorial **AS**sociation **MAP**ping

A repository with the source code of the CASMAP package (Application Note submission).

**Contents:**
+ [Downloading the repository](#downloading-the-repository)
+ [Installation and compilation](#installation-and-compilation)
+ [Examples](#examples)

## Downloading the repository

Clone the repository with the git command
```
git clone https://github.com/BorgwardtLab/CASMAP.git
```
The location of the downloaded repository will be referred to as *root folder*.

Alternatively, download the repository as a ZIP file and decompress in your local computer.


## Installation and compilation

**Firstly, all dependencies must be installed (R, Python, various packages - see detailed instructions below).** After that, installing CASMAP is simple.

### Step 1: Obtaining the C++ compiler

The current version uses the GCC 6 compiler (`gcc-6`) in the Makefiles.
**For Linux (Ubuntu)**, GCC 6 can be installed [using the following commands](https://askubuntu.com/questions/746369/how-can-i-install-and-use-gcc-6-on-xenial/746480#746480):
```
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-6 g++-6
```

**For Mac OSX**, install GCC 6 with the following commands:
```
brew update
brew install gcc@6
```


### Step 2: Compiling the code base and R/Python packages

The commands below apply to both **Linux** and **Mac OSX**.

Change directory to the *root folder*. There, simply run:

#### Compilation of core functions

```
make
```

#### Compilation of the R wrapper

From the *root folder* change directory to the subdirectory `r_wrapper`:

```
cd r_wrapper
```

Run the command:

```
make package
```

**Troubleshooting**: If the previous command does not work, then run:

```
Rscript -e "devtools::build()"
```

which will then create the R package .tar.gz file.


#### Compilation of the Python wrapper

From the *root folder* change directory to the subdirectory `cython_wrapper`:

```
cd cython_wrapper
```

Run the command:

```
make
```

**Troubleshooting**: If the previous command fails with `make: *** [wrapper] Error 1`, then run:

```
export CC=gcc-6
make
```

which will make sure that the compilation is done with a newer version of the `gcc` compiler.


#### Uninstall

To remove the compiled files (perhaps before recompiling), in the *root folder* run:

```
make clean
```

### Step 3: Installing the R package

Follow the steps described in the [README.md](r_wrapper/README.md) file under the subdirectory `r_wrapper`.

### Step 4: Installing the Python package

Follow the steps described in the [README.md](cython_wrapper/README.md) file under the subdirectory `cython_wrapper`.


## Troubleshooting

On **Ubuntu 16.04**, make sure the following packages are installed:

 * libxml2-dev (needed for R's devtools)
 * curl (needed for R's devtools)
 * libcurl4-openssl-dev (needed for R's devtools)
 * libssl-dev (needed for R's devtools)
 * ipython (for convenience)
 * cython (for the python wrapper)
 * python-pip (installs setuptools, which is needed for the python wrapper)

These can be installed with the following commands.

```
sudo apt-get install libxml2-dev
sudo apt-get install curl
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libssl-dev
sudo apt install python-pip
sudo apt install ipython
sudo apt install cython
```


### The R wrapper

In order to have use the R wrapper, [the R language](https://cran.r-project.org/) needs to be installed. Moreover, a few packages need to be installed, in particular `devtools`. In R:


```
install.packages("devtools", dependencies=T)
```

Note that this will take a long time! If a few dependencies are not successfully installed, just make sure that the following works in R:


```
library(devtools)
```
Follow the steps described in the README.md file under the subdirectory `r_wrapper`.


### The Python wrapper

Similarly, for the Python wrapper to work, we need Python to be installed. We also need cython and pip to be installed. One way to do this is to follow the instructions [on this page](http://pip.readthedocs.io/en/stable/installing/).

### Installing sphinx

If one uses

```
make clean
```

in the cython folder, the message that `sphinx-build` is missing might appear. To install it in Ubuntu:

```
sudo apt-get install python-sphinx
```

To install it in OSX:

```
sudo pip install --upgrade pip
pip install sphinx --user
```

Follow the steps described in the README.md file under the subdirectory `cython_wrapper`.

## Examples

The repository includes [examples](examples) that illustrate how to use CASMAP to carry out region-based association studies and higher-order epistasis search.

The folder [examples/data](examples/data) provides a real-world dataset from the plant model organism *A. thaliana*, downloaded from the [easyGWAS](https://easygwas.ethz.ch/) online resource.

+ [examples/data/region_based/avrB](examples/data/region_based/avrB) contains: 1) `X.dat`, the genotypes of 87 samples measured at 214,032 homozygous SNPs; 2) Y.dat, the binary phenotype *avrB* for each of the 87 samples, 3) `C.dat`, a categorical covariate to account for population structure; and 4) `plink.map`, a list with the locations of all SNPs in the genome of *A. thaliana*.
+ [examples/data/higher_order_epistasis/avrB](examples/data/higher_order_epistasis/avrB) is analogous, but file `X.dat` in this folder consists of a subset of 650 SNPs in Chromosome 1.

The folder [examples/code](examples/code) contains Python/R Jupyter notebooks which show how CASMAP can be used to analyze these datasets.

+  `run_region_based.ipynb` (Python) and `run_region_based.Rmd` (R) describe the region-based association study example.
+ `run_higher_order_epistasis.ipynb` (Python) and `run_higher_order_epistasis.Rmd` (R) describe the higher-order epistasis search example.

The notebooks are heavily commented, aiming to provide a comprehensive, guided explanation of how to use CASMAP. For additional details about these examples, we refer the user to the Supplementary Material of the Application Note.

Additional [examples](https://www.ethz.ch/content/dam/ethz/special-interest/bsse/borgwardt-lab/Projects/ISMB18-tutorial/casmap_tutorial.zip), as well as [slides](https://www.ethz.ch/content/dam/ethz/special-interest/bsse/borgwardt-lab/Projects/ISMB18-tutorial/module1.pdf) describing the theory behind the methods included in CASMAP, can be found as part of the material of the following [tutorial](https://www.bsse.ethz.ch/mlcb/education/tutorial-ismb18.html), which took place in ISMB 2018.
