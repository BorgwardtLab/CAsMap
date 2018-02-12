# README

Contents of this folder:

* `libsigpatsearch.doxyfile` - Doxygen configuration file;
* `Makefile` - to facilitate build and cleanup.

## Prerequisites

Install `doxygen` and `dot` command line tools from, respectively Doxygen
(>= 1.8) and GraphViz software packages. You can skip the latter, and prevent
graphs generation by adjusting the Doxygen configuration file.

## Build

To build docs run in the current folder:

    $ make

This will built in the default configuration the HTML docs in the `html/`
subfolder of this folder. You can alter the output target (incl. format) in the
configuration file.

## Read

Open main target file (by default `html/index.html`).

## Cleanup

To remove previously built docs run in the current folder:

    $ make clean

This will remove target folder (by default `html/`) with all of its contents.
