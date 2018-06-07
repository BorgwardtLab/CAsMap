# CASMAP
Combinatorial Association MAPping

A repository with the source code of the CASMAP package (Application Note submission).

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

Troubleshooting: If the previous command does not work, then run:

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

#### Uninstall

To remove the compiled files (perhaps before recompiling), in the *root folder* run:

```
make clean
```

### Step 3: Installing the R package

Follow the steps described in the README.md file under the subdirectory `r_wrapper`.

### Step 4: Installing the Python package

Follow the steps described in the README.md file under the subdirectory `cython_wrapper`.


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

Note that before one could simply use

```
import sigpatsearch as s
```

but now for Ubuntu it seems it is necessary to use 

```
import sigpatsearch.wrapper as s
```

in order to access the ```createSigPatSearch``` command from the new wrapper/constructor, as in:

```
sig_fais1 = s.createSigPatSearch(method='fais')
```


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
