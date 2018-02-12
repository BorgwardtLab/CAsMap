# README

## Installation

### Prerequisites

You need a C++ compiler that sufficiently supports a C++11 standard, i.e.
GCC >= 5, or, recommended, LLVM >= 7; cf.

    $ g++ --version

You need `python` and `pip` installed on your computer.

On a Linux OS, you may need to separately install Python's header files, e.g.
with `apt` (Ubuntu, or Debian):

    $ sudo apt-get install python-dev

or with `yum` (CentOS, Red Hat, Fedora):

    $ sudo yum install python-devel

To install requirements for the wrappers execute from this folder:

    $ make install_requirements

which should install all tools needed to build Cython wrappers and to test them.

*Note*: for pip installation of packages Makefile forces by default `--user`
flag (w/ no path). You can override that behaviour by setting or passing
directly to `make` command `PIP_INSTALL_FLAGS="..."` variable.

### Build, install and test

To build the wrappers run from this folder:

    $ make all

To install the wrappers from the folder run:

    $ make install

Optionally, to generate separate HTML docs for the installed packae run:

   $ make doc_html

Note that the output folder given is relative to the `doc/` folder.

To test the installed package run:

    $ make test

**Beware**: tests check output which depends on using `rand`/`srand` from
`stdlib` (WY methods), so despite setting same seed you might still get
_slightly_ different results between different versions of `stdlib` (different
OSs in particular). Checks should always pass for non-WY methods and with
C libraries provided by Apple LLVM version 8.0.0.

#### Troubleshooting

Make sure that you have successfully installed prerequisites, in particular that the
`pytest` and `sphinx` binaries are visible; try:

    $ which py.test
    $ which sphinx-build

If these return empty strings, instead of paths to the installed binaries, then modify
your `PATH` environmental variable to include paths of the local installation binaries.

If your packages end up in a wrong Python's version folder (e.g. ver. 3 instead of the
default ver. 2), then use an explicit pip binary matching you default Python version,
e.g.:

    $ make install_requirements all install PIP=pip2

### Build, package and install from the package

Alternatively, you can build, package and install from the package by running:

    $ make install_pkg

### Cleanup

To remove build and packaging files run:

    $ make clean

To uninstall the package from your system run:

    $ make uninstall
