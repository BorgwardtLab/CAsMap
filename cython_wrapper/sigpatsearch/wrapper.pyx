from wrapper cimport *
from libcpp cimport bool

from warnings import warn
from collections import namedtuple

#for wrapper
import numbers
import math
import warnings


MemoryUsage = namedtuple("MemoryUsage", ["current", "peak"])



IntervalWithP = namedtuple("IntervalWithP", ["start", "end", "pvalue"])

Region = namedtuple("Region", "start end")

_SummaryBase = namedtuple("Summary", "testability_threshold significance_level corrected_signficance_level")
_SummaryInt = namedtuple("Summary", ('int_processed', 'int_testable',) + _SummaryBase._fields)
_SummaryFais = namedtuple('Summary', _SummaryInt._fields + ('testability_region',))


_ResultInt = namedtuple("Result", "summary sig_int sig_int_clustered")
class ResultInt(_ResultInt):

    def __str__(self):
        res = _ResultInt(
            self.summary,
            "<{} significant intervals>".format(len(self.sig_int)),
            "<{} clustered intervals>".format(len(self.sig_int_clustered)),
        )
        return str(res)


cdef class _SignificantFeaturesSearch:
    cdef object alpha
    cdef object lmax
    cdef bool _results_available
    cdef bool _file_loaded

    def __cinit__(self, set_defaults=True):
        """Create a new instance.

        :param set_defaults: Set default alpha and lmax values (default: True)
        """
        self.alpha = self.lmax = None
        self._results_available = self._file_loaded = False
        if set_defaults:
            self.set_lmax(0)
            self.set_alpha(0.05)

    def _check_if_alpha_value_is_allowed(self, double x):
        assert x >= 0 and x <= 1, "you need to set alpha to a value between 0 and 1"

    def _set_results_unavailable(self):
       self._results_available = False
    def _set_results_available(self):
       self._results_available = True

    def set_alpha(self, double alpha):
        """Set the significance threshold (FWER).

        :param alpha: Significance threshold value
        """
        self._check_if_alpha_value_is_allowed(alpha)
        self.alpha = alpha
        self._set_results_unavailable()

    def get_alpha(self):
        """Get the significance threshold (FWER).
        """
        return self.alpha

    def _check_if_lmax_value_is_allowed(self, x):
        assert isinstance(x, (int, long,)) and x >= 0, "you need to set lmax to a non-negative integer value"

    def set_lmax(self, lmax):
        """Set the maximum number of SNPs to consider.

        :param lmax:  Max. number of SNPs (0 means all in the file)
        """
        self._check_if_lmax_value_is_allowed(lmax)
        self.lmax = lmax
        self._set_results_unavailable()

    def get_lmax(self):
        """Set the maximum number of SNPs to consider.
        """
        return self.lmax

    def read_eth_files(self, str x_file, str y_file):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        """
        self._check_if_read_is_allowed()
        try:
            self._readETHFiles(x_file, y_file)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def read_plink_files(self, str base_file_path):
        """Read the PLINK format files.

        :param: base_file_path: base name path string (w/o extension) for all PLINK files
        """
        self._check_if_read_is_allowed()
        try:
            self._readPlinkFiles(base_file_path)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def _check_if_files_are_loaded(self):
        assert self._file_loaded, "you have to load a data file first"

    def _check_if_read_is_allowed(self):
        pass

    def write_eth_files(self, str x_file, str y_file):
        """Write ETH format files.

        :param x_file: Data file path string
        :param y_file: Labels file path string
        """
        self._check_if_write_is_allowed()
        self._writeETHFiles(x_file, y_file)

    def _check_if_write_is_allowed(self):
        self._check_if_files_are_loaded()

    def execute(self):
        """Execute the search.
        """
        self._check_if_execute_is_allowed()
        self._set_results_unavailable()
        self._execute(self.alpha, self.lmax)
        self._set_results_available()

    def _check_if_alpha_is_set(self):
        assert self.alpha is not None, "you have to call set_alpha first"

    def _check_if_lmax_is_set(self):
        assert self.lmax is not None, "you have to call set_lmax first"

    def _check_if_execute_is_allowed(self):
        self._check_if_files_are_loaded()
        self._check_if_alpha_is_set()
        self._check_if_lmax_is_set()

    def _check_if_results_available(self):
        assert self._results_available, "you need to call the execute method first"

    def write_summary(self, str path):
        """Write search summary to a file.

        :param path: Write file path string
        """
        self._check_if_results_available()
        self._write_summary(path)

    def write_profile(self, str path):
        """Write profiling summary to a file.

        :param path: Write file path string
        """
        self._check_if_results_available()
        self._write_profile(path)

    def get_libmem(self):
        """Return current and peak memory (RSS) usage of the underlying library.

        :return `MemoryUsage` namedtuple.
        """
        return MemoryUsage(
            self._get_curr_memory(),
            self._get_peak_memory(),
        )


    #NOTE NEW: print function
    def __repr__(self):
        myclasstype = str(type(self).__name__)
        mytype = "unknown"
        if myclasstype=="_SignificantIntervalSearchExact":
            mytype = "FAIS"
        if myclasstype=="_SignificantIntervalSearchFastCmh":
            mytype = "FastCMH"
        if myclasstype=="_SignificantItemsetSearchFacs":
            mytype = "FACS"

        #extract alpha and lmxax
        myalpha = str(self.alpha)
        mylmax = str(self.lmax)

        #output message to be returned
        message = mytype + " object with:\n" + \
                  " * alpha = " + myalpha + "\n" + \
                  " * lmax = " + mylmax
        return(message)


cdef class _SignificantIntervalSearch(_SignificantFeaturesSearch):
    def write_filtered_intervals(self, str path):
        """Write filtered intervals (most significant from each cluster) to a file.

        :param path: Write file path string
        """
        self._check_if_results_available()
        self._write_filtered_intervals(path)

    def write_pvals_testable_intervals(self, str path):
        """Write p-values for all testable intervals found to a file.

        :param path: Write file path string
        """
        self._check_if_results_available()
        self._write_pvals_testable_intervals(path)

    def write_pvals_significant_intervals(self, str path):
        """Write p-values for all significant intervals found to a file.

        :param path: Write file path string
        """
        self._check_if_results_available()
        self._write_pvals_significant_intervals(path)

    cdef _intervals(self, vector[Interval] ivec):
        cdef list intervals = []
        cdef Interval iv
        cdef int i
        for i in range(ivec.size()):
            iv = ivec.at(i)
            intervals.append(IntervalWithP(iv.getStart(), iv.getEnd(), iv.getPvalue()))
        return intervals

    def get_significant_intervals(self):
        """Returns found significant intervals.

        :return Vector of intervals
        """
        self._check_if_results_available()
        return self._get_significant_intervals()

    def get_filtered_intervals(self):
        """Returns most significant intervals from each cluster of overlapping
        significant intervals.

        :return Vector of intervals
        """
        self._check_if_results_available()
        return self._get_filtered_intervals()

    def get_summary(self):
        """Return search summary.

        :return `Summary` namedtuple.
        """
        self._check_if_results_available()
        return self._get_summary_int()

    def get_result(self):
        """Return search results.

        :return `Result` namedtuple.
        """
        self._check_if_results_available()

        return ResultInt(self.get_summary(), self.get_significant_intervals(),
                      self.get_filtered_intervals())

cdef class _SignificantIntervalSearchFais(_SignificantIntervalSearch):
    def get_summary(self):
        """Return search summary.

        :return `Summary` namedtuple.
        """
        self._check_if_results_available()
        return self._get_summary_fais()

cdef class _SignificantIntervalSearchExact(_SignificantIntervalSearchFais):
    """Exact fast significant interval search

    Class for exact significant intervals search with Taron correction for
    bounding intermediate FWERs.
    """
    cdef SignificantIntervalSearchExact inst

    def __cinit__(self, **kwargs):
        self.inst = SignificantIntervalSearchExact()

    def _readETHFiles(self, x_file, y_file):
        self.inst.readETHFiles(x_file, y_file)

    def _readPlinkFiles(self, base_file_path):
        self.inst.readPlinkFiles(base_file_path)

    def _writeETHFiles(self, x_file, y_file):
        self.inst.writeETHFiles(x_file, y_file)

    def _execute(self, alpha, lmax):
        self.inst.execute(alpha, lmax)

    def _write_summary(self, path):
        self.inst.getSummary().writeToFile(path)

    def _write_profile(self, path):
        self.inst.getProfiler().writeToFile(path)

    def _write_filtered_intervals(self, path):
        self.inst.getFilteredIntervals().writeToFile(path)

    def _write_pvals_testable_intervals(self, path):
        self.inst.getPValsTestableInts().writeToFile(path)

    def _write_pvals_significant_intervals(self, path):
        self.inst.getPValsSigInts().writeToFile(path)

    def _get_significant_intervals(self):
        cdef vector[Interval] ivec = self.inst.getSignificantIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_filtered_intervals(self):
        cdef vector[Interval] ivec = self.inst.getFilteredIntervals().getSigInts()
        return self._intervals(ivec)


    def _get_summary_fais(self):
        cdef SummaryFais summary = self.inst.getSummary()
        region = Region((summary.getSl1(), summary.getSu1()), (summary.getSl2(), summary.getSu2()))

        return _SummaryFais(
            summary.getNumFeaturesProcessed(), summary.getm(),
            summary.getDelta(), summary.getAlpha(), summary.getDelta_opt(),
            region
        )

    def _get_curr_memory(self):
        return self.inst.getProfiler().getCurrMemory()

    def _get_peak_memory(self):
        return self.inst.getProfiler().getPeakMemory()

cdef class _SignificantIntervalSearchChi(_SignificantIntervalSearchFais):

    cdef SignificantIntervalSearchChi inst

    def __cinit__(self, **kwargs):
        self.inst = SignificantIntervalSearchChi()

    def _readETHFiles(self, x_file, y_file):
        self.inst.readETHFiles(x_file, y_file)

    def _readPlinkFiles(self, base_file_path):
        self.inst.readPlinkFiles(base_file_path)

    def _writeETHFiles(self, x_file, y_file):
        self.inst.writeETHFiles(x_file, y_file)

    def _execute(self, alpha, lmax):
        self.inst.execute(alpha, lmax)

    def _write_summary(self, path):
        self.inst.getSummary().writeToFile(path)

    def _write_profile(self, path):
        self.inst.getProfiler().writeToFile(path)

    def _write_filtered_intervals(self, path):
        self.inst.getFilteredIntervals().writeToFile(path)

    def _write_pvals_testable_intervals(self, path):
        self.inst.getPValsTestableInts().writeToFile(path)

    def _write_pvals_significant_intervals(self, path):
        self.inst.getPValsSigInts().writeToFile(path)

    def _get_significant_intervals(self):
        cdef vector[Interval] ivec = self.inst.getSignificantIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_filtered_intervals(self):
        cdef vector[Interval] ivec = self.inst.getFilteredIntervals().getSigInts()
        return self._intervals(ivec)


    def _get_summary_fais(self):
        cdef SummaryFais summary = self.inst.getSummary()
        region = Region((summary.getSl1(), summary.getSu1()), (summary.getSl2(), summary.getSu2()))

        return _SummaryFais(
            summary.getNumFeaturesProcessed(), summary.getm(),
            summary.getDelta(), summary.getAlpha(), summary.getDelta_opt(),
            region
        )

    def _get_curr_memory(self):
        return self.inst.getProfiler().getCurrMemory()

    def _get_peak_memory(self):
        return self.inst.getProfiler().getPeakMemory()


cdef class _SignificantIntervalSearchWy(_SignificantIntervalSearchFais):
    cdef object n_perm

    def __cinit__(self, set_defaults=True, **kwargs):
        """Create a new instance.

        :param set_defaults: Set default alpha, lmax and n_perm values (default: True)
        """
        self.n_perm = None
        if (set_defaults):
            self.set_n_perm(50)

    def _check_if_n_perm_value_is_allowed(self, x):
        assert isinstance(x, (int, long,)) and x > 0, "you need to set n_perm to a positive integer value"

    def set_n_perm(self, n_perm):
        """Set number of permutations of the observations used for estimation of intermediate FWERs.

        :param n_perm: Number of permutations (positive integer)
        """
        self._check_if_n_perm_value_is_allowed(n_perm)
        self._set_n_perm(n_perm)
        self.n_perm = n_perm
        self._set_results_unavailable()

    def get_n_perm(self, n_perm):
        """Get number of permutations of the observations used for estimation of intermediate FWERs.
        """
        return self.n_perm

    def set_seed(self, unsigned seed):
        """Set random seed (for reproducibility).

        :param seed: Random number generator seed (non-negative integer)
        """
        self._set_seed(seed)
        self._set_results_unavailable()

cdef class _SignificantIntervalSearchWyChi(_SignificantIntervalSearchWy):

    cdef SignificantIntervalSearchWyChi inst

    def __cinit__(self, **kwargs):
        self.inst = SignificantIntervalSearchWyChi()

    def _readETHFiles(self, x_file, y_file):
        self.inst.readETHFiles(x_file, y_file)

    def _readPlinkFiles(self, base_file_path):
        self.inst.readPlinkFiles(base_file_path)

    def _writeETHFiles(self, x_file, y_file):
        self.inst.writeETHFiles(x_file, y_file)

    def _execute(self, alpha, lmax):
        self.inst.execute(alpha, lmax)

    def _write_summary(self, path):
        self.inst.getSummary().writeToFile(path)

    def _write_profile(self, path):
        self.inst.getProfiler().writeToFile(path)

    def _write_filtered_intervals(self, path):
        self.inst.getFilteredIntervals().writeToFile(path)

    def _write_pvals_testable_intervals(self, path):
        self.inst.getPValsTestableInts().writeToFile(path)

    def _write_pvals_significant_intervals(self, path):
        self.inst.getPValsSigInts().writeToFile(path)

    def _get_significant_intervals(self):
        cdef vector[Interval] ivec = self.inst.getSignificantIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_filtered_intervals(self):
        cdef vector[Interval] ivec = self.inst.getFilteredIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_summary_fais(self):
        cdef SummaryFais summary = self.inst.getSummary()
        region = Region((summary.getSl1(), summary.getSu1()), (summary.getSl2(), summary.getSu2()))

        return _SummaryFais(
            summary.getNumFeaturesProcessed(), summary.getm(),
            summary.getDelta(), summary.getAlpha(), summary.getDelta_opt(),
            region
        )

    def _get_curr_memory(self):
        return self.inst.getProfiler().getCurrMemory()

    def _get_peak_memory(self):
        return self.inst.getProfiler().getPeakMemory()

    def _set_n_perm(self, n_perm):
        self.inst.setNPerm(n_perm)

    def _set_seed(self, seed):
        self.inst.setSeed(seed)


cdef class _SignificantIntervalSearchWyExact(_SignificantIntervalSearchWy):

    cdef SignificantIntervalSearchWyExact inst

    def __cinit__(self, **kwargs):
        self.inst = SignificantIntervalSearchWyExact()

    def _readETHFiles(self, x_file, y_file):
        self.inst.readETHFiles(x_file, y_file)

    def _readPlinkFiles(self, base_file_path):
        self.inst.readPlinkFiles(base_file_path)

    def _writeETHFiles(self, x_file, y_file):
        self.inst.writeETHFiles(x_file, y_file)

    def _execute(self, alpha, lmax):
        self.inst.execute(alpha, lmax)

    def _write_summary(self, path):
        self.inst.getSummary().writeToFile(path)

    def _write_profile(self, path):
        self.inst.getProfiler().writeToFile(path)

    def _write_filtered_intervals(self, path):
        self.inst.getFilteredIntervals().writeToFile(path)

    def _write_pvals_testable_intervals(self, path):
        self.inst.getPValsTestableInts().writeToFile(path)

    def _write_pvals_significant_intervals(self, path):
        self.inst.getPValsSigInts().writeToFile(path)

    def _get_significant_intervals(self):
        cdef vector[Interval] ivec = self.inst.getSignificantIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_filtered_intervals(self):
        cdef vector[Interval] ivec = self.inst.getFilteredIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_summary_fais(self):
        cdef SummaryFais summary = self.inst.getSummary()
        region = Region((summary.getSl1(), summary.getSu1()), (summary.getSl2(), summary.getSu2()))

        return _SummaryFais(
            summary.getNumFeaturesProcessed(), summary.getm(),
            summary.getDelta(), summary.getAlpha(), summary.getDelta_opt(),
            region
        )

    def _get_curr_memory(self):
        return self.inst.getProfiler().getCurrMemory()

    def _get_peak_memory(self):
        return self.inst.getProfiler().getPeakMemory()

    def _set_n_perm(self, n_perm):
        self.inst.setNPerm(n_perm)

    def _set_seed(self, seed):
        self.inst.setSeed(seed)



# Note: no multiple inheritance in Cython, so this part is duplicated for isets
cdef class _SignificantIntervalSearchWithCovariates(_SignificantIntervalSearch):
    cdef bool _cov_loaded

    def __cinit__(self, **kwargs):
        self._cov_loaded = False

    def _check_if_execute_is_allowed(self):
        self._check_if_files_are_loaded()
        self._check_if_covariates_are_loaded()
        self._check_if_alpha_is_set()
        self._check_if_lmax_is_set()

    def _check_if_covariates_are_loaded(self):
        if self._cov_loaded:
            return
        warn("assuming one covariate for all observations; to change covariates call the update_covariates_file method first")

    def read_eth_files(self, str x_file, str y_file, object cov_file = None):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readETHFilesWithCovariates(x_file, y_file, cov_file)
                self._cov_loaded = True
            else:
                self._readETHFiles(x_file, y_file)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def read_plink_files(self, str base_file_path, object cov_file = None):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readPlinkFilesWithCovariates(base_file_path, cov_file)
                self._cov_loaded = True
            else:
                self._readPlinkFiles(base_file_path)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def update_covariates_file(self, str cov_file):
        """Update covariates file (in the same format as ETH labels file). If not
        read, by default, a single covariate for all observations is assumed.

        :param cov_file: Covariates file path
        """
        self._check_if_files_are_loaded()
        try:
            self._update_covariates_file(cov_file)
        except RuntimeError as e:
            raise IOError(e.message)
        self._cov_loaded = True
        self._set_results_unavailable()

    def write_covariates_file(self, str cov_file):
        """Write covariates files.

        :param cov_file: Covariates file path
        """
        self._write_covariates_file(cov_file)



cdef class _SignificantIntervalSearchFastCmh(_SignificantIntervalSearchWithCovariates):

    cdef SignificantIntervalSearchFastCmh inst

    def __cinit__(self, **kwargs):
        self.inst = SignificantIntervalSearchFastCmh()

    def _readETHFiles(self, x_file, y_file):
        self.inst.readETHFiles(x_file, y_file)

    def _readPlinkFiles(self, base_file_path):
        self.inst.readPlinkFiles(base_file_path)

    def _readETHFilesWithCovariates(self, x_file, y_file, cov_file):
        self.inst.readETHFilesWithCovariates(x_file, y_file, cov_file)

    def _readPlinkFilesWithCovariates(self, base_file_path, cov_file):
        self.inst.readPlinkFilesWithCovariates(base_file_path, cov_file)

    def _writeETHFiles(self, x_file, y_file):
        self.inst.writeETHFiles(x_file, y_file)

    def _execute(self, alpha, lmax):
        self.inst.execute(alpha, lmax)

    def _write_summary(self, path):
        self.inst.getSummary().writeToFile(path)

    def _write_profile(self, path):
        self.inst.getProfiler().writeToFile(path)

    def _write_filtered_intervals(self, path):
        self.inst.getFilteredIntervals().writeToFile(path)

    def _write_pvals_testable_intervals(self, path):
        self.inst.getPValsTestableInts().writeToFile(path)

    def _write_pvals_significant_intervals(self, path):
        self.inst.getPValsSigInts().writeToFile(path)

    def _get_significant_intervals(self):
        cdef vector[Interval] ivec = self.inst.getSignificantIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_filtered_intervals(self):
        cdef vector[Interval] ivec = self.inst.getFilteredIntervals().getSigInts()
        return self._intervals(ivec)

    def _get_curr_memory(self):
        return self.inst.getProfiler().getCurrMemory()

    def _get_peak_memory(self):
        return self.inst.getProfiler().getPeakMemory()

    def _get_summary_int(self):
        cdef SummaryInt summary = self.inst.getSummary()

        return _SummaryInt(
            summary.getNumFeaturesProcessed(), summary.getm(),
            summary.getDelta(), summary.getAlpha(), summary.getDelta_opt()
        )

    def _update_covariates_file(self, cov_file):
        self.inst.readCovariatesFile(cov_file)

    def _write_covariates_file(self, cov_file):
        self.inst.writeCovariatesFile(cov_file)



ItemsetWithP = namedtuple("ItemsetWithP", ["itemset", "pvalue"])

_SummaryFacs = namedtuple("Summary", ('iset_processed', 'iset_closed_processed', 'iset_testable',) + _SummaryBase._fields)

_ResultIset = namedtuple("Result", "summary sig_iset")
class ResultIset(_ResultIset):

    def __str__(self):
        res = _ResultIset(
            self.summary,
            "<{} significant itemsets>".format(len(self.sig_iset)),
        )
        return str(res)



cdef class _SignificantItemsetSearch(_SignificantFeaturesSearch):

    def write_pvals_testable_itemsets(self, str path):
        """Write p-values for all testable itemsets found to a file.

        :param path: Write file path string
        """
        self._check_if_results_available()
        self._write_pvals_testable_itemsets(path)

    def write_pvals_significant_itemsets(self, str path):
        """Write p-values for all significant itemsets found to a file.

        :param path: Write file path string
        """
        self._check_if_results_available()
        self._write_pvals_significant_itemsets(path)

    cdef _itemsets(self, ItemsetSet iset):
        cdef vector[vector[longint]] ivec = iset.getItemsetsVector();
        cdef vector[double] pvals = iset.getPValueVector();

        cdef list itemsets = []
        cdef int i
        for i in range(ivec.size()):
            itemsets.append(ItemsetWithP(ivec.at(i), pvals.at(i)))
        return itemsets

    def get_significant_itemsets(self):
        """Returns found significant itemsets.

        :return Vector of itemsets
        """
        self._check_if_results_available()
        return self._get_significant_itemsets()

    def get_summary(self):
        """Return search summary.

        :return `Summary` namedtuple.
        """
        self._check_if_results_available()
        return self._get_summary_facs()
        # Rem: should be _get_summary_iset(), but we don't have usage for it,
        #      hence FACS override directly here

    def get_result(self):
        """Return search results.

        :return `Result` namedtuple.
        """
        self._check_if_results_available()
        return ResultIset(self.get_summary(), self.get_significant_itemsets())

# Note: no multiple inheritance in Cython, so this part is duplicated for ints
cdef class _SignificantItemsetSearchWithCovariates(_SignificantItemsetSearch):
    cdef bool _cov_loaded

    def __cinit__(self, **kwargs):
        self._cov_loaded = False

    def _check_if_execute_is_allowed(self):
        self._check_if_files_are_loaded()
        self._check_if_covariates_are_loaded()
        self._check_if_alpha_is_set()
        self._check_if_lmax_is_set()

    def _check_if_covariates_are_loaded(self):
        if self._cov_loaded:
            return
        warn("assuming one covariate for all observations; to change covariates call the update_covariates_file method first")

    def read_eth_files(self, str x_file, str y_file, object cov_file = None):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readETHFilesWithCovariates(x_file, y_file, cov_file)
                self._cov_loaded = True
            else:
                self._readETHFiles(x_file, y_file)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def read_plink_files(self, str base_file_path, object cov_file = None):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readPlinkFilesWithCovariates(base_file_path, cov_file)
                self._cov_loaded = True
            else:
                self._readPlinkFiles(base_file_path)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def update_covariates_file(self, str cov_file):
        """Update covariates file (in the same format as ETH labels file). If not
        read, by default, a single covariate for all observations is assumed.

        :param cov_file: Covariates file path
        """
        self._check_if_files_are_loaded()
        try:
            self._update_covariates_file(cov_file)
        except RuntimeError as e:
            raise IOError(e.message)
        self._cov_loaded = True
        self._set_results_unavailable()

    def write_covariates_file(self, str cov_file):
        """Write covariates files.

        :param cov_file: Covariates file path
        """
        self._write_covariates_file(cov_file)



cdef class _SignificantItemsetSearchFacs(_SignificantItemsetSearchWithCovariates):

    cdef SignificantItemsetSearchFacs inst

    def __cinit__(self, **kwargs):
        self.inst = SignificantItemsetSearchFacs()


    def _readETHFiles(self, x_file, y_file):
        self.inst.readETHFiles(x_file, y_file)

    def _readPlinkFiles(self, base_file_path):
        self.inst.readPlinkFiles(base_file_path)

    def _readETHFilesWithCovariates(self, x_file, y_file, cov_file):
        self.inst.readETHFilesWithCovariates(x_file, y_file, cov_file)

    def _readPlinkFilesWithCovariates(self, base_file_path, cov_file):
        self.inst.readPlinkFilesWithCovariates(base_file_path, cov_file)

    def _writeETHFiles(self, x_file, y_file):
        self.inst.writeETHFiles(x_file, y_file)

    def _execute(self, alpha, lmax):
        self.inst.execute(alpha, lmax)

    def _write_summary(self, path):
        self.inst.getSummary().writeToFile(path)

    def _write_profile(self, path):
        self.inst.getProfiler().writeToFile(path)

    def _write_pvals_testable_itemsets(self, path):
        self.inst.getPValsTestableIsets().writeToFile(path)

    def _write_pvals_significant_itemsets(self, path):
        self.inst.getPValsSigIsets().writeToFile(path)

    def _get_significant_itemsets(self):
        return self._itemsets(self.inst.getPValsSigIsets())

    def _get_curr_memory(self):
        return self.inst.getProfiler().getCurrMemory()

    def _get_peak_memory(self):
        return self.inst.getProfiler().getPeakMemory()

    def _get_summary_facs(self):
        cdef SummaryFacs summary = self.inst.getSummary()

        return _SummaryFacs(
            summary.getNumFeaturesProcessed(),
            summary.getNumItemsetsClosedProcessed(), summary.getm(),
            summary.getDelta(), summary.getAlpha(), summary.getDelta_opt()
        )

    def _update_covariates_file(self, cov_file):
        self.inst.readCovariatesFile(cov_file)

    def _write_covariates_file(self, cov_file):
        self.inst.writeCovariatesFile(cov_file)



def isInOpenInterval(x, lower=0, upper=1):
    '''Checks if a value is numeric and strictly between two other values.

       Args:
            x (numeric): Value to be checked. Needs to be numeric.
            lower (numeric): Lower bound. Default value is '0'.
            upper (numeric): Upper bound. Default value is '1'.

       Returns:
            If numeric, and  strictly greater than 'lower' and 
            strictly smaller than 'upper', then return 'True'. 
            Else return 'False'.
    '''
    inInterval = True
    if isinstance(x, numbers.Number):
        if  math.isnan(x) | math.isinf(x):
            inInterval = False
        else:
            if (x <= lower) | (x >= upper): 
                inInterval = False
    else:
        #not numeric
        inInterval = False

    return(inInterval)


def checkIsBoolean(var, name):
    '''Checks if a variable is a boolean, if not throws error
       otherwise returns boolean

       Args: 
            var (boolean?): The variable to be checked (if boolean).
            name (string): The name of the variable to appear in any 
                           error message.

       Returns: 
            If not boolean (or 'None'), throws error.
            If 'None', returns 'False'.
            Otherwise, return value of boolean.

    '''
    if var is not None:
        #now check boolean or not
        if var in [True, False]:
            return(var)
        else:
            message = "Error: " + str(name) + " is not a boolean."
            raise ValueError(message)
    else:
        #is None, return False
        return(False)



def createSigPatSearch(*args, 
                       method='', 
                       use_intervals=None, 
                       use_combinations=None,
                       use_covariate=False,
                       alpha=0.05,
                       max_length=0,
                       **kwargs):
    '''Creates an object to run one of the following methods:
        * FAIS
        * FastCMH
        * FACS
       
       There are two ways to initilise the object, either using the
       'method' parameter or using the 'use_intervals' boolean flag.
       More details below.

       Args:
        * 'method': one of 'fais', 'fastcmh' or 'facs' (upper/lower case 
          does not matter). Default is empty string ''.

        * 'use_intervals': boolean for setting whether intervals are used 
                            (i.e. True for FAIS and FastCMH), 
                            or not (i.e. False for FACS)

        * 'use_combinations': boolean for setting whether combinations are 
                              used (i.e. True for FACS) or not  
                              (i.e. False for FAIS and FACS). This is the 
                              opposite of use_intervals, but note that if
                              there is a conflict, in that both use_intervals
                              and use_combinations are set to the same value,
                              a method with intervals will be used
                              (FAIS or FastCMH).

        * 'use_covariate': boolean for setting whether a covariate is used
                            (i.e. True for FastCMH or FACS) or not 
                            (i.e. False for FAIS). Default value is False.

        * 'alpha': A numerical value strictly between 0 and 1 which specifies
                   the significance threshold. Default value is 0.05.

        * 'max_length': The length of the largest interval/combination that
                        is searched for, i.e. if max_length=5, then only
                        intervals/combinations of length 1, 2, 3, 4 and 5 will
                        be considered. Setting this value to 0 means ALL
                        possible lengths are considered. Default value is 0.

        Examples to follow...

    '''

    #check maxlength
    # Python does not seem to have a nice way to check for finite
    # numbers, when argument could be string, nan or inf
    if isinstance(max_length, numbers.Number):
        if  math.isnan(max_length) | math.isinf(max_length):
            message = "'max_length' needs to be either 0 or a positive integer."
            raise ValueError(message)
        else:
            max_length = int(max_length)
            if max_length < 0:
                max_length = 0
    else:
        message = "'max_length' needs to be either 0 or a positive integer."
        raise ValueError(message)


    #check alpha
    if not isInOpenInterval(alpha):
        message ="'alpha' needs to be a value strictly in interval (0, 1)."
        raise ValueError(message)


    #check alpha
    if not isInOpenInterval(alpha):
        message ="'alpha' needs to be a value strictly in interval (0, 1)."
        raise ValueError(message)


    #TODO: add support for LAMP (use_combinations==True and use_covariate=False)
    #TODO: Add examples
    #TODO: check default case
    #TODO: Add check for unrecognized arguments; for example if "blah=False" is
    #      passed to the function
    if (method=='') & (use_intervals is None) & (use_combinations is None):
        message = "Need to set at least one of 'method',"
        message += " 'use_intervals' or 'use_combinations'."
        raise ValueError(message)

    if (method != ''):
        #convert to lower case (first convert to string)
        method = str(method).lower()


        #check if one of 'fais', 'fastcmh' or 'facs'
        correctName = False

        #fais
        if method == 'fais':
            correctName = True
            use_intervals = True
            use_combinations = False
            use_covariate = False
            sig = _SignificantIntervalSearchExact()
            sig.set_alpha(alpha) 
            sig.set_lmax(max_length)
            return(sig)


        #fastcmh
        if method == 'fastcmh':
            correctName = True
            use_intervals = True
            use_combinations = False
            use_covariate = True
            sig = _SignificantIntervalSearchFastCmh()
            sig.set_alpha(alpha) 
            sig.set_lmax(max_length)
            return(sig)


        #facs
        if method == 'facs':
            correctName = True
            use_intervals = False
            use_combinations = True
            use_covariate = True
            sig = _SignificantItemsetSearchFacs()
            sig.set_alpha(alpha) 
            sig.set_lmax(max_length)
            return(sig)

        if correctName == False:
            message = "'method' parameter needs to be one of " + \
                      "'fais', 'fastcmh' or 'facs'. Otherwise, " + \
                      "set the 'use_intervals' or 'use_combinations' " + \
                      "flags."
            raise ValueError(message)

    #end of if (method != '')

    #now we look at the case that the name as NOT been set, but the flags
    #use_intervals or use_combinations have been set.

    #First, check how many 'None's
    none_intervals = (use_intervals is None)
    none_combinations = (use_combinations is None)
    trueCount = none_intervals + none_combinations
    if trueCount==2:
        message = "Need to set 'method', or at least one of the " + \
                  "boolean flags 'use_intervals' and " + \
                  "'use_combinations'."
        raise ValueError(message)

    #NOW check flags are all booleans
    #and convert 'None' values to 'False' values
    use_intervals = checkIsBoolean(use_intervals, "use_intervals")
    use_combinations = checkIsBoolean(use_combinations, "use_combinations")
    use_covariate = checkIsBoolean(use_covariate, "use_covariate")


    # At this point, either use_intervals or use_combinations must be NOT 'None'
    # 'None' values receive lower priority
    if none_intervals:
        use_intervals = not use_combinations

    if none_combinations:
        use_combinations = not use_intervals


    #do fais
    if (use_intervals | (not use_combinations)) & (not use_covariate):
        sig = _SignificantIntervalSearchExact()
        sig.set_alpha(alpha) 
        sig.set_lmax(max_length)
        return(sig)
    
    #do fastcmh
    if (use_intervals | (not use_combinations)) & use_covariate:
        sig = _SignificantIntervalSearchFastCmh()
        sig.set_alpha(alpha) 
        sig.set_lmax(max_length)
        return(sig)
    
    #do facs/lamp
    if use_combinations | (not use_intervals):
        sig = _SignificantItemsetSearchFacs()
        sig.set_alpha(alpha) 
        sig.set_lmax(max_length)
        return(sig)

    #should never reach this point
    #if here, return error
    message = "Error: no correct flags were used."
    raise ValueError(message)

    #end of createSigPatSearch

