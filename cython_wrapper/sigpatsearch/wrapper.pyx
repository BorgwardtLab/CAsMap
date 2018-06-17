from wrapper cimport *
from libcpp cimport bool

from warnings import warn
from collections import namedtuple

#for wrapper
import numbers
import math
import warnings


MemoryUsage = namedtuple("MemoryUsage", ["current", "peak"])



IntervalWithP = namedtuple("IntervalWithP", ["start", "end", "score", "odds_ratio", "pvalue"])

Region = namedtuple("Region", "start end")

_SummaryBase = namedtuple("Summary", "testability_threshold target_fwer corrected_significance_threshold")
_SummaryInt = namedtuple("Summary", ('n_int_processed', 'n_int_testable',) + _SummaryBase._fields)
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

    def read_eth_files(self, str x_file, str y_file, str encoding):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        """
        self._check_if_read_is_allowed()
        try:
            self._readETHFiles(x_file, y_file, encoding)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def read_plink_files(self, str base_file_path, str encoding):
        """Read the PLINK format files.

        :param: base_file_path: base name path string (w/o extension) for all PLINK files
        """
        self._check_if_read_is_allowed()
        try:
            self._readPlinkFiles(base_file_path, encoding)
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
        if myclasstype=="_SignificantIntervalSearchExact": # TODO: What about Chi2?
            mytype = "FAIS"
        if myclasstype=="_SignificantIntervalSearchFastCmh":
            mytype = "FastCMH"
        if myclasstype=="_SignificantItemsetSearchFacs":
            mytype = "FACS"

        #extract alpha and lmxax
        myalpha = str(self.alpha)
        mylmax = str(self.lmax)

        #output message to be returned:
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
            intervals.append(IntervalWithP(iv.getStart(), iv.getEnd(), iv.getScore(), iv.getOddsRatio(), iv.getPvalue()))
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

    def _readETHFiles(self, x_file, y_file, encoding):
        self.inst.readETHFiles(x_file, y_file, encoding)

    def _readPlinkFiles(self, base_file_path, encoding):
        self.inst.readPlinkFiles(base_file_path, encoding)

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

    def _readETHFiles(self, x_file, y_file, encoding):
        self.inst.readETHFiles(x_file, y_file, encoding)

    def _readPlinkFiles(self, base_file_path, encoding):
        self.inst.readPlinkFiles(base_file_path, encoding)

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
        #warn("assuming one covariate for all observations; to change covariates call the update_covariates_file method first")

    def read_eth_files(self, str x_file, str y_file, object cov_file = None, str encoding = "dominant"):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        :param encoding: < "dominant" | "recessive" >. Default: "dominant"
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readETHFilesWithCovariates(x_file, y_file, cov_file, encoding)
                self._cov_loaded = True
            else:
                self._readETHFiles(x_file, y_file, encoding)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def read_plink_files(self, str base_file_path, object cov_file = None, str encoding = "dominant"):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        :param encoding: < "dominant" | "recessive" >. Default: "dominant"
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readPlinkFilesWithCovariates(base_file_path, cov_file, encoding)
                self._cov_loaded = True
            else:
                self._readPlinkFiles(base_file_path, encoding)
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

    def _readETHFiles(self, x_file, y_file, encoding):
        self.inst.readETHFiles(x_file, y_file, encoding)

    def _readPlinkFiles(self, base_file_path, encoding):
        self.inst.readPlinkFiles(base_file_path, encoding)

    def _readETHFilesWithCovariates(self, x_file, y_file, cov_file, encoding):
        self.inst.readETHFilesWithCovariates(x_file, y_file, cov_file, False, encoding)

    def _readPlinkFilesWithCovariates(self, base_file_path, cov_file, encoding):
        self.inst.readPlinkFilesWithCovariates(base_file_path, cov_file, True, encoding)

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



ItemsetWithP = namedtuple("ItemsetWithP", ["itemset", "score", "odds_ratio", "pvalue"])

_SummaryFacs = namedtuple("Summary", ('n_iset_processed', 'n_iset_closed_processed', 'n_iset_testable',) + _SummaryBase._fields)

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

    cdef _itemsets(self, ItemsetSetWithOddsRatio iset):
        cdef vector[vector[longint]] ivec = iset.getItemsetsVector();
        cdef vector[double] scores = iset.getScoreVector();
        cdef vector[double] odds_ratios = iset.getOddsRatioVector();
        cdef vector[double] pvals = iset.getPValueVector();

        cdef list itemsets = []
        cdef int i
        for i in range(ivec.size()):
            itemsets.append(ItemsetWithP(ivec.at(i), scores.at(i), odds_ratios.at(i), pvals.at(i)))
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
        #warn("assuming one covariate for all observations; to change covariates call the update_covariates_file method first")

    def read_eth_files(self, str x_file, str y_file, object cov_file = None, str encoding = "dominant"):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        :param encoding: < "dominant" | "recessive" >. Default: "dominant"
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readETHFilesWithCovariates(x_file, y_file, cov_file, encoding)
                self._cov_loaded = True
            else:
                self._readETHFiles(x_file, y_file, encoding)
        except RuntimeError as e:
            raise IOError(e.message)
        self._set_results_unavailable()
        self._file_loaded = True

    def read_plink_files(self, str base_file_path, object cov_file = None, str encoding = "dominant"):
        """Read the ETH format files.

        :param x_file: Data file path
        :param y_file: Labels file path
        :param cov_file: Covariates file path (default: None, means single covariate for all observations)
        :param encoding: < "dominant" | "recessive" >. Default: "dominant"
        """
        self._check_if_read_is_allowed()
        try:
            if cov_file != None:
                self._readPlinkFilesWithCovariates(base_file_path, cov_file, encoding)
                self._cov_loaded = True
            else:
                self._readPlinkFiles(base_file_path, encoding)
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


    def _readETHFiles(self, x_file, y_file, encoding):
        self.inst.readETHFiles(x_file, y_file, encoding)

    def _readPlinkFiles(self, base_file_path, encoding):
        self.inst.readPlinkFiles(base_file_path, encoding)

    def _readETHFilesWithCovariates(self, x_file, y_file, cov_file, encoding):
        print('caca')
        self.inst.readETHFilesWithCovariates(x_file, y_file, cov_file, False, encoding)

    def _readPlinkFilesWithCovariates(self, base_file_path, cov_file, encoding):
        self.inst.readPlinkFilesWithCovariates(base_file_path, cov_file, 1, encoding)

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


class CASMAP(object):

    _ALLOWABLE_MODES = ('regionGWAS', 'higherOrderEpistasis')
    _ALLOWABLE_ENCODINGS = ('dominant', 'recessive')

    def __init__(self, mode, alpha=0.05, max_comb_size=0):
        self.setMode(mode)
        self.setTargetFWER(alpha)
        self.setMaxCombinationSize(max_comb_size)

        self._core = None
        self._use_covariates = None

    def getMode(self):
        return self._mode

    def getTargetFWER(self):
        return self._alpha

    def getMaxCombinationSize(self):
        return self._max_comb_size

    def isInitialized(self):
        return self._core is not None

    def _checkInitialized(self):
        if self._core is None:
            raise ValueError('Object not initialized or hyperparameters changed since last execution. Please call method readFiles prior to execute.')

    def setMode(self, mode):
        # check mode
        if mode not in CASMAP._ALLOWABLE_MODES:
            raise ValueError("Currently implemented modes: < " + " | ".join(CASMAP._ALLOWABLE_MODES) + " >.")
        self._mode = mode

        # Delete previous "core" object (if any)
        self._core = None

    def setTargetFWER(self, alpha=0.05):
        #check alpha
        if not isInOpenInterval(alpha):
            raise ValueError("Target FWER 'alpha' needs to be a value strictly between 0 and 1.")
        self._alpha = alpha

        # Delete previous "core" object (if any)
        self._core = None

    def setMaxCombinationSize(self, max_comb_size=0):
        #check maxlength
        # Python does not seem to have a nice way to check for finite
        # numbers, when argument could be string, nan or inf
        if isinstance(max_comb_size, numbers.Number):
            if  math.isnan(max_comb_size) | math.isinf(max_comb_size):
                raise ValueError("Maximum combination size 'max_comb_size' needs to be either 0 (unlimited) or a positive integer.")
            else:
                self._max_comb_size = int(max_comb_size)
                if self._mode == 'higherOrderEpistasis' and self._max_comb_size > 0:
                    print("The current implementation of higher-order epistasis analyses does not support a limited maximum number of interacting variants. The analysis will be carried out for an unlimited order.")
                    self._max_comb_size = 0
                if self._max_comb_size < 0:
                    self._max_comb_size = 0
        else:
            raise ValueError("Maximum combination size 'max_comb_size' needs to be either 0 (unlimited) or a positive integer.")

        # Delete previous "core" object (if any)
        self._core = None

    def _createCore(self):
        if self._use_covariates is not None:
            # Instantiate object of the appropriate depending on options
            if self._mode == 'regionGWAS' and self._use_covariates is False:
                self._core = _SignificantIntervalSearchChi()
            elif self._mode == 'regionGWAS' and self._use_covariates is True:
                self._core = _SignificantIntervalSearchFastCmh()
            elif self._mode == 'higherOrderEpistasis':
                self._core = _SignificantItemsetSearchFacs()

            # Set parameters of the object
            self._core.set_alpha(self._alpha)
            self._core.set_lmax(self._max_comb_size)
        else:
            self._core = None

    def readFiles(self, genotype_file=None, phenotype_file=None, plink_file_root=None, covariate_file=None, encoding="dominant"):
        # Check whether user decided to use tab-separated text files (binary_format) or PLINK formatted files (plink_format)
        binary_format = genotype_file is not None and phenotype_file is not None
        plink_format = plink_file_root is not None
        # At least one of the two must be two, otherwise raise an error
        if not (binary_format or plink_format):
            raise ValueError('Either plink_file_root or genotype_file and phenotype_file must be specified as arguments.')
        # Check that encoding type is correct
        if encoding not in CASMAP._ALLOWABLE_ENCODINGS:
            raise ValueError("Currently implemented encodings: < " + " | ".join(CASMAP._ALLOWABLE_ENCODINGS) + " >.")

        # If an additional covariates file was specified, set the object into "CMH mode"
        self._use_covariates = covariate_file is not None

        # Create appropriate "core" object
        self._createCore()

        # Give preference to plink_format over binary_format if, by any reason, a user decides to mess around and
        # specify both
        if plink_format:
            if self._use_covariates:
                self._core.read_plink_files(plink_file_root, covariate_file, encoding)
            else:
                self._core.read_plink_files(plink_file_root, encoding)
        elif binary_format:
            if self._use_covariates:
                self._core.read_eth_files(genotype_file, phenotype_file, covariate_file, encoding)
            else:
                self._core.read_eth_files(genotype_file, phenotype_file, encoding)
        else:  # this branch should not be reachable due to the check above...
            raise ValueError('Either plink_file_root or genotype_file and phenotype_file must be specified as arguments.')

    def execute(self):
        self._checkInitialized()
        self._core.execute()

    def writeSummary(self, path):
        self._checkInitialized()
        self._core.write_summary(path)

    def writeProfile(self, path):
        self._checkInitialized()
        self._core.write_profile(path)

    def writeSignificantRegions(self, path):
        if self._mode != 'regionGWAS':
            raise ValueError('Method writeSignificantRegions only available for region-based GWAS analyses.')
        self._checkInitialized()
        self._core.write_pvals_significant_intervals(path)

    def writeSignificantClusterRepresentatives(self, path):
        if self._mode != 'regionGWAS':
            raise ValueError('Method writeSignificantClusterRepresentatives only available for region-based GWAS analyses.')
        self._checkInitialized()
        self._core.write_filtered_intervals(path)

    def writeSignificantInteractions(self, path):
        if self._mode != 'higherOrderEpistasis':
            raise ValueError('Method writeSignificantInteractions only available for higher-order epistasis analyses.')
        self._checkInitialized()
        self._core.write_pvals_significant_itemsets(path)

    def getSummary(self):
        self._checkInitialized()
        return self._core.get_summary()

    def getSignificantRegions(self):
        if self._mode != 'regionGWAS':
            raise ValueError('Method getSignificantRegions only available for region-based GWAS analyses.')
        self._checkInitialized()
        return self._core.get_significant_intervals()

    def getSignificantClusterRepresentatives(self):
        if self._mode != 'regionGWAS':
            raise ValueError('Method getSignificantClusterRepresentatives only available for region-based GWAS analyses.')
        self._checkInitialized()
        return self._core.get_filtered_intervals()

    def getSignificantInteractions(self):
        if self._mode != 'higherOrderEpistasis':
            raise ValueError('Method getSignificantInteractions only available for higher-order epistasis analyses.')
        self._checkInitialized()
        return self._core.get_significant_itemsets()

    def __repr__(self):
        message = "CASMAP object with:\n" + \
                  " * Mode = {}".format(self._mode) + "\n" + \
                  " * Target FWER = {}".format(self._alpha) + "\n" + \
                  " * Maximum combination size = {}".format(self._max_comb_size) + "\n"
        if self._core is not None:
            message += " * Input files read\n" + \
                       " * Covariate = {}".format(self._use_covariates)
        else:
            message += " * No input files read\n"

        return message
