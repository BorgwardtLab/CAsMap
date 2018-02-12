Significant Pattern Search C++ library {#mainpage}
======================================

This is developers documentation for the Significant Pattern Search C++
library. You can navigate here the documentation generated from the source code.
Additionally, below you will find technical documentation describing particular
higher-level aspects of the library.

> Attention: with GNU C++ Compiler at least version 5 is required for some
> of the C++11 standard library features which were not fully implemented in
> lower versions of the compiler.

[TOC]

Adding or modifying an input format {#addinput}
===================================

Data (SignificantPattern::Genotype class), as well as labels and covariates
(SignificantPattern::Phenotype class) are all stored in memory as
`unsigned char` array or vector, respectively. The common abstraction for types
of inputs is an abstract SignificantPattern::ArrayFile class.

> Note: `unsigned char` gives a implementation-independent minimal memory
> footprint (`sizeof(char) == 1`; cf.
> [`bool`](http://en.cppreference.com/w/cpp/language/types#Boolean_type) for
> binary data), but, on the other hand, it imposes also an limitation of max.
> 256 values, specifically for non-binary covariates.

What happens on read?
---------------------

Data and labels files are read simultaneously. To invoke read user has to call
either of the API read methods for the supported data formats:
SignificantPattern::SignificantFeaturesSearch::readETHFiles(), or
SignificantPattern::SignificantFeaturesSearch::readPlinkFiles(). These
functions call the protected
SignificantPattern::SignificantFeaturesSearch::readFiles()
function which (indirectly) is responsible for:

* checking and reading into the temporary buffer the small labels file,
* checking and reading into the object's memory the large data file,
* memoization of the labels file if no error occurred while reading the data
  file, or rolling back to the memory state from before the read,
* timing I/O,
* delegating actual checks and read to the data and labels classes.

Eventually, one of the format-specific methods  `readFORMATFile()` on the labels
(SignificantPattern::Phenotype) or data (SignificantPattern::Genotype) class
is called, where `FORMAT=ETH,PlnikFam,PlinkCov,PlinkRaw`. Each of these methods
does two runs through the file calling two private methods:

1. first, `checkFORMATFile()`, which does not read full file content into the
   memory, but reads chuncks of file to check for errors in the file format, and
   to read input parameters such as data size or presence of a header;
2. then, `parseFORMATFile()`, which does actual read file contents into the
   memory.

For checking a size consistency, the labels size if passed to the data read
method. Additionally, because the PLINK format does not assume same ordering of
the labels and data, mapping from the line identifier to the row numbers read
from the labels file (SignificantPattern::Phenotype::getPlinkLineForFIDAndIID)
is passed to and used while reading a PLINK covariates or data file to check
inputs consistency and to set memory contents in a row-aligned fashion (like
in the ETH format).

How about covariates?
---------------------

Currently, for methods supporting covariates (subclasses of
SignificantPattern::SignificantFeaturesSearchWithCovariates class), the
covariates file:

* is supported in ETH file format, or [PLINK .cov short
  format](https://www.cog-genomics.org/plink2/formats#cov) (3 columns)
* can be read together with data and labels files, or after these files have
  been read first (cf.
  SignificantPattern::SignificantFeaturesSearchWithCovariates::readETHFilesWithCovariates,
  SignificantPattern::SignificantFeaturesSearchWithCovariates::readPlinkFilesWithCovariates,
  or SignificantPattern::SignificantFeaturesSearchWithCovariates::readCovariatesFile),

So how do I add or modify a read method?
----------------------------------------

To modify an existing read method for a file format `FORMAT`, modify the
corresponding the `checkFORMATFile()` and `parseFORMATFile()` methods on
SignificantPattern::Genotype or SignificantPattern::Phenotype class, and, if
required the `readFORMATFile()` method, which calls them both.

To add a read method for a new format, except for adding the above-mentioned
methods on SignificantPattern::Genotype and SignificantPattern::Phenotype
classes, you should also make the read available in the API class
SignificantPattern::SignificantIntervalSearch by adding a corresponding
`readFORMATFiles()` method and extending the protected
SignificantPattern::SignificantFeaturesSearch::readFiles() method which the
former should always call.

And what about write?
---------------------

Write is done in ETH format with calls to the API methods:
SignificantPattern::SignificantFeaturesSearch::writeETHFiles(), or
SignificantPattern::SignificantFeaturesSearchWithCovariates::writeCovariatesFile()
and
SignificantPattern::SignificantFeaturesSearchWithCovariates::writeETHFilesWithCovariates().
All of these methods delegate the write to SignificantPattern::Phenotype and
SignificantPattern::Genotype classes, which, to that end, implement
a pure (abstract) virtual method
SignificantPattern::ArrayFile::writeFileStream().

Adding a new search method or extending a current one {#addmethod}
=====================================================

Executing pattern search can be done only once the input files have been read.
The non-virtual API method is
SignificantPattern::SignificantFeaturesSearch::execute(),
which itself calls multiple virtual methods, which are then overridden,
depending on the search method. Please see the inheritance diagram for the
SignificantPattern::SignificantIntervalSearch class for an overview of
available search methods.

What happens on search execution?
---------------------------------

SignificantPattern::SignificantFeaturesSearch::execute() method is responsible
for:

* main search flow: initializing a search, invoking the search for patterns,
  post-processing results of the search, and cleaning up after the search,
* timing subsequent parts of the search, and marking the memory usage,
* cleanup in case of an error.

### 1. Initialization

Initialization related to a single search execution is done in the
SignificantPattern::SignificantFeaturesSearch::execute_init() method, which
breaks down into the steps as described next.

#### 1.1 Pre-search cleanup and zero-initialization

> Beware: this part is subject to adjustments in the near future.

First, pre-search cleanup is required for when results from the previous search
are still in the memory. Second, zero-initialization is a concept of setting
manually-managed pointers to `0` to be able to recognize if they have been
initialized, but here also calling default (empty) constructors for
automatically-managed dynamic data structures.

Done in two calls:

1. a virtual method `execute_destructor()`, which first invokes private method
   `execute_destructor_METHOD()` which performs cleanup operations local to the
   current class of search, and then makes a call to the super class.
2. a virtual method `execute_constructor()`, which first makes a call to the
   super class, and then invokes a private method `execute_constructor_METHOD()`
   which performs zero-initialization local to the current class of search.

> Note: this approach mimics a default destructor, and a default constructor
> calls, as they would have been done automatically when the search-specific
> parts of the memory would be separated (from inputs read/write) into its own
> class.

#### 1.2 General variables initialization

For instance, allocation of memory for a temporary parent layer data buffer.
Done directly in
SignificantPattern::SignificantFeaturesSearch::execute_init().

#### 1.3 Algorithm-specific initialization

Allocation and initialization of temporary data specific to the search
class representing a method (an algorithm), such as a  frequency counters for
intervals,  frequency counters for intervals, log-gamma cache, permutations, or
pre-computed values per each coviariate class. Done in a virtual method
`algorithm_init()`.

### 2. Search
Search is done in two rounds of search: first to compute a corrected
significance threshold (p-value), second to find the actual patterns, given the
corrected threshold.

#### 2.1. Computing corrected significant threshold

Done in a virtual method `compute_corrected_significance_threshold()`.

In most of the methods there is additional structure to this method, which is to
invoke `process_first_layer_threshold()` virtual method, and then
`process_intervals_threshold()` virtual method.

> Note: currently FACS method is setup greedily, and omits the first round (see
> SignificantPattern::SignificantItemsetSearchFacs::compute_corrected_significance_threshold())

#### 2.2 Finding significant patterns

Done in a virtual method `find_significant_features()`.

In most of the currently implemented methods there is additional structure to
this method, which is to invoke `process_first_layer_pvalues()` virtual
method, and then `process_intervals_pvalues()` virtual method.

### 3. Algorithm-specific cleanup and post-processing

Algorithm-specific cleanup and results is done in the
SignificantPattern::SignificantFeaturesSearch::execute_end() method, which
breaks down into the steps as described next.

#### 3.1 Post-search algorithm-specific cleanup

Cleanup of the data set up before in the corresponding `algorithm_init()`
virtual method.

#### 3.2 Setting up output features

Done in `process_significant_features()` virtual method

#### 3.3 Setting up summary

Done directly in SignificantPattern::SignificantFeaturesSearch::execute_end().

So how do I add a new search method?
------------------------------------

Depending on how similar a new method is to the already implemented methods, you
have to pick a place in the search class hierarchy for it, and implement the
virtual methods mentioned in the search execution description above (according
to their purpose). If necessary, you might need to re-factor the current code,
to introduce a more general code structure for the significant patterns search.

Testing and debugging {#testdebug}
=====================

Sanity check regression tests
-----------------------------

To test the library you have to build the provisional executables (built when
running `make` from the top-level folder), and then run one or more Bash
regression test; from the top-level folder:

    $ test/regression_test.sh # short ETH-format based tests for all methods
    $ test/regression_test_plink.sh # long PLINK-format based tests for all methods

To runs tests selectively adjust `STEM_*` variables in the above-mentioned Bash
scripts.

> Note: on Mac OS (Darwin kernel) you have to "code sign" debugger to be able to
> use it.

In case of errors set `ECHO_CMD=true` in a test script to print the commands
that are actually executed.

Debugging
---------

Before debugging, re-build all sources using the `DEBUG=1` make variable, i.e.:

    $ make clean all DEBUG=1

This will not only leave non-optimized symbols in the source code, but also will
print a lot of debugging output. When re-running a failing executable you should
be able to roughly identify the location of an error by just reading the output.

To debug with a debugger use your favorite: by default `gdb` or, if compiled
with LLVM, `lldb`; then run, e.g.:

    $ lldb -- executables/significant_itemset_search_facs -eth data/facs/test/tictactoe_data.txt data/facs/test/tictactoe_labels.txt data/facs/test/tictactoe_covariates.txt 0.05 0 test/output/output_iset_facs

> Note: on Mac OS (Darwin kernel) you have to "code sign" debugger to be able to
> use it.

For debugger commands see on-line docs e.g. http://lldb.llvm.org/lldb-gdb.html .

Extending wrappers {#wrappers}
==================

Cython wrapper
--------------

Approach: encapsulate library (concrete) classes. This allows to hook
pre-/post- lib’s calls (early checks, process errors or results). Encapsulation
is done using Cython classes (in `.pyx` source files), which import directly
from C++ (in `.pxd` header files).

In `cython_wrapper` folder:

1. symlink new library sources in `sigpatsearch/`
2. edit `sigpatsearch/wrapper.pxd`, then `sigpatsearch/wrapper.pyx` files
   analogously to existing methods (with docstrings)
3. add regression and unit tests in `tests/test_wrappers.py` analogously to
   existing tests
   * see [`pytest-regtest`](https://pypi.python.org/pypi/pytest-regtest) for
     instruction on how to regress tests

R wrapper
---------

Approach: encapsulate library methods as separate functions. Encapsulation is
done using Rcpp and RefClasses for `obj$method()` syntax.

> Note: R wrapper should eventually directly use of C++ classes with Rcpp
> modules, analogously to the Cython wrapper. Multi-inheritance may be an issue
> though for Rcpp.

> Note: RefClasses implement pass-by-reference (mutable) vs. R’s usual
> copy-on-modify semantics and implements message-passing OO (methods to call
> are passed as a part of the call), so methods belong to classes, not functions
> (enables `obj$method()` use).

In `r_wrapper/` folder:

1. symlink missing headers and sources in `src/` and headers in `inst/include/`
2. edit Rcpp wrappers in `src/sigpatsearch.cpp`, then R encapsulation classes in
   `R/wrapper.r`,  and docs in `man-roxygen/*.R` files, analogusly to existing
   methods
3. add unit tests in `test/testthat.R`
