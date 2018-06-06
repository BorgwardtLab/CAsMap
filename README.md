# CAsMap

A repository with the software for the Application Note submission.

### Data files

 **Note:** The large data files have been excluded here - please do not add large data files to the repo, at least for the moment.

Note that the data folder can be downloaded from the original [SIS repo on GitLab](https://sissource.ethz.ch/sispub/significant_interval_search_rewrite)

Specifically, the data folder is [here](https://sissource.ethz.ch/sispub/significant_interval_search_rewrite/tree/master/data)

## Installation/compilation

If all of the dependencies are installed (R, Python, various packages - see below), then installation is simple. In the root folder, simply run


```
make
```

### R wrapper compilation

Move to the folder

```
cd r_wrapper
```

and use the following command:


```
make package
```

If that does not work in Linux, then run:

```
Rscript -e "devtools::build()"
```

which will then create the R package .tar.gz file.


### Python wrapper compilation

Similarly, in the folder

```
cd cython_wrapper
```

run

```
make
```

### Testing

For both the R and Python wrappers, run:

```
make test
```

in the respective folders.




### Uninstall

To remove the compiled files (perhaps before recompiling), in the root directory run:

```
make clean
```



## Troubleshooting


On Ubuntu 16.04, make sure the following packages are installed:

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

### Using gcc6

The current version uses the gcc6 compiler in the makefile, in an attempt to remove minor differences between Linux and OSX. It is not necessary, but I think there were a few minor numerical differences using gcc5. Anyway, gcc6 can be installed [using the following commands](https://askubuntu.com/questions/746369/how-can-i-install-and-use-gcc-6-on-xenial/746480#746480):

```
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-6 g++-6
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

