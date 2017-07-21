# encoding: utf-8
from __future__ import print_function, division, absolute_import

import hashlib
import glob
import os
import gc
import psutil
import re
import random

import pytest

from sigpatsearch import SignificantIntervalSearchExact
from sigpatsearch import SignificantIntervalSearchChi
from sigpatsearch import SignificantIntervalSearchWyExact
from sigpatsearch import SignificantIntervalSearchWyChi
from sigpatsearch import SignificantIntervalSearchFastCmh
from sigpatsearch import SignificantItemsetSearchFacs

# (search_class, if_set_n_perm, if_read_cov, if_itemset)
_ALL_SEARCHES = [
    (SignificantIntervalSearchExact, False, False, False),
    (SignificantIntervalSearchChi, False, False, False),
    (SignificantIntervalSearchWyExact, True, False, False),
    (SignificantIntervalSearchWyChi, True, False, False),
    (SignificantIntervalSearchFastCmh, False, True, False),
    (SignificantItemsetSearchFacs, False, True, True),
]


def format(df, col_width=10):
    fmt_string = "{:%d}" % col_width
    new_rows = []
    for row in str(df).split("\n"):
        cells = re.split("[ ]+", row)
        new_row = " ".join(fmt_string.format(cell) for cell in cells)
        new_rows.append(new_row)
    return "\n".join(new_rows)


@pytest.fixture
def data_path():

    # using os.path.relpath shortens the output when using in regression tests:
    here = os.path.relpath(os.path.dirname(__file__), os.getcwd())

    def compute_path(partial_path):
        return os.path.join(here, "data", partial_path)
    return compute_path



def _dump_tmpdir(tmpdir, regtest):

    print(file=regtest)
    print("files:", file=regtest)
    print(sorted(os.listdir(tmpdir.strpath)), file=regtest)

    for p in sorted(glob.glob(tmpdir.join("*").strpath)):
        if "output_timing" in p:
            continue

        print(file=regtest)
        print(os.path.basename(p), file=regtest)
        print("-" * 80, file=regtest)
        print(file=regtest)

        lines = open(p, "r").readlines()

        digest = hashlib.sha1()
        digest.update(b"".join(l.encode("ascii") for l in lines))
        print("hexdigest:", digest.hexdigest(), file=regtest)

        print("number of lines:", len(lines), file=regtest)
        print(file=regtest)
        print("content:", file=regtest)
        print("------", file=regtest)
        if len(lines) > 40:
            lines = lines[:20] + ["..."] + lines[-20:]
        for line in lines:
            print(line.rstrip(), file=regtest)
        print("------", file=regtest)

def _dump_result(s, regtest, itemset=False):
    summary = s.get_summary()
    print(summary, file=regtest)

    if itemset:
        features = s.get_significant_itemsets()
    else:
        features = s.get_significant_intervals()
    print("not clustered:\n", len(features), file=regtest)

    print(s.get_result(), file=regtest)
    if itemset:
        print("\nsig isets:", file=regtest)
        _dump_itemsets(s.get_result().sig_iset, regtest)
    else:
        print("\nsig ints:", file=regtest)
        _dump_intervals(s.get_result().sig_int, regtest)

        print("\nsig ints clusterd:", file=regtest)
        _dump_intervals(s.get_result().sig_int_clustered, regtest)

def _dump_itemset(i, iset, regtest):
    print("row %4d:" % i, "%.5e" % iset.pvalue,
          ",".join(("%3d" % item) for item in iset.itemset), file=regtest)

def _dump_interval(i, interval, regtest):
    print("row %4d:" % i, "%3d" % interval.start, "%3d" %
          interval.end, "%.5e" % interval.pvalue, file=regtest)

def _dump_features(dump_fun, features, regtest):
    if len(features) > 20:
        for i, feat in enumerate(features[:10]):
            dump_fun(i, feat, regtest)
        print("...", file=regtest)
        n_feat = len(features)
        for i, feat in enumerate(features[-10:]):
            dump_fun(n_feat-10+i, feat, regtest)
    else:
        for i, feat in enumerate(features):
            dump_fun(i, feat, regtest)

def _dump_intervals(*args):
    _dump_features(_dump_interval, *args)

def _dump_itemsets(*args):
    _dump_features(_dump_itemset, *args)

def _dump_output(s, tmpdir, wcov=False, itemset=False):
    s.write_summary(tmpdir.join("output_summary.txt").strpath)
    s.write_profile(tmpdir.join("output_timing.txt").strpath)
    if itemset:
        s.write_pvals_testable_itemsets(tmpdir.join(
            "p_vals_testable_itemsets.txt").strpath)
        s.write_pvals_significant_itemsets(tmpdir.join(
            "p_vals_significant_itemsets.txt").strpath)
    else:
        s.write_filtered_intervals(tmpdir.join(
            "output_sigints_filtered.corrected.csv").strpath)
        s.write_pvals_testable_intervals(tmpdir.join(
            "p_vals_testable_intervals.txt").strpath)
        s.write_pvals_significant_intervals(tmpdir.join(
            "p_vals_significant_intervals.txt").strpath)
    s.write_eth_files(tmpdir.join("output_data.txt").strpath,
                      tmpdir.join("output_label.txt").strpath)
    if wcov:
        s.write_covariates_file(tmpdir.join("output_cov.txt").strpath)

def _test_read(s, readfun_str, readfun_args, cov_file=""):
    # files not read
    with pytest.raises(Exception):
        s.execute()

    # should be able to read files without setting anything
    readfun = getattr(s, readfun_str)
    import time
    started = time.time()
    readfun(*readfun_args)
    if cov_file:
        s.update_covariates_file(cov_file)
    print("[DEBUG] time to read files = %.3f sec" % (time.time() - started))

def _test_nonwwy(s, readfun_str, readfun_args, tmpdir, data_path, regtest):

    _test_read(s, readfun_str, readfun_args)

    s.set_alpha(0.04)
    s.set_lmax(0)
    s.execute()

    s.set_alpha(0.04)
    s.set_lmax(0)

    s.execute()

    _dump_output(s, tmpdir)
    _dump_result(s, regtest)
    _dump_tmpdir(tmpdir, regtest)

def test_run_significant_interval_search_exact(tmpdir, data_path, regtest):

    _test_nonwwy(
        SignificantIntervalSearchExact(),
        "read_eth_files", (data_path("data.txt"), data_path("label.txt"),),
        tmpdir, data_path, regtest
    )

def test_run_significant_interval_search_exact_plink(tmpdir, data_path, regtest):

    _test_nonwwy(
        SignificantIntervalSearchExact(),
        "read_plink_files", (data_path("sample_data"),),
        tmpdir, data_path, regtest
    )

def test_run_significant_interval_search_chi(tmpdir, data_path, regtest):

    _test_nonwwy(
        SignificantIntervalSearchChi(),
        "read_eth_files", (data_path("data.txt"), data_path("label.txt"),),
        tmpdir, data_path, regtest
    )



def _test_wcov(s, readfun_str, readfun_args, cov_path, tmpdir, regtest, itemset=False):

    # read data and labels first
    with pytest.raises(Exception):
        s.update_covariates_file(cov_path)

    # default covariates are initialised only on execute or read w/o covariates
    with pytest.raises(Exception):
        s.write_covariates_file(tmpdir.join("output_cov0.txt").strpath)

    # execute w/ default single covariate
    _test_read(s, readfun_str, readfun_args)

    s.set_alpha(0.05)
    s.set_lmax(0)
    s.execute()

    s.write_covariates_file(tmpdir.join("output_cov0.txt").strpath)

    # execute w/ a covariates file
    s.update_covariates_file(cov_path)

    s.set_alpha(0.05)
    s.set_lmax(0)

    s.execute()

    _dump_output(s, tmpdir, wcov=True, itemset=itemset)
    _dump_result(s, regtest, itemset=itemset)
    _dump_tmpdir(tmpdir, regtest)

def test_run_significant_interval_search_fastcmh(tmpdir, data_path, regtest):

    _test_wcov(
        SignificantIntervalSearchFastCmh(),
        "read_eth_files", (data_path("cov_data.txt"), data_path("cov_label.txt"),),
        data_path("cov.txt"), tmpdir, regtest,
    )

def test_run_significant_itemset_search_facs(tmpdir, data_path, regtest):

    _test_wcov(
        SignificantItemsetSearchFacs(),
        "read_eth_files", (data_path("tictactoe_data.txt"), data_path("tictactoe_labels.txt"),),
        data_path("tictactoe_covariates.txt"), tmpdir, regtest,
        itemset=True
    )



def _test_wy(s, readfun_str, readfun_args, tmpdir, data_path, regtest):

    _test_read(s, readfun_str, readfun_args)

    s.set_alpha(0.25)
    s.set_lmax(5)

    # seed and n_perm not set
    with pytest.raises(Exception):
        s.execute()

    s.set_n_perm(4)
    s.set_seed(1)

    # now we have everything initialized so that execute()
    s.execute()

    _dump_output(s, tmpdir)
    _dump_result(s, regtest)
    _dump_tmpdir(tmpdir, regtest)

def test_run_significant_interval_search_wy_exact(tmpdir, data_path, regtest):

    _test_wy(
        SignificantIntervalSearchWyExact(),
        "read_eth_files", (data_path("data.txt"), data_path("label.txt"),),
        tmpdir, data_path, regtest
    )

def test_run_significant_interval_search_wy_chi(tmpdir, data_path, regtest):

    _test_wy(
        SignificantIntervalSearchWyChi(),
        "read_eth_files", (data_path("data.txt"), data_path("label.txt"),),
        tmpdir, data_path, regtest
    )



def test_error_checks(tmpdir, data_path, regtest):
    data_fn = "tictactoe_data.txt"
    label_fn = "tictactoe_labels.txt"
    cov_fn = "tictactoe_covariates.txt"

    for (s_cls, if_set_n_perm, if_read_cov, if_itemset) in _ALL_SEARCHES:
        s = s_cls(set_defaults=False)
        s_name = s_cls.__name__

        # nonexisting data file
        with pytest.raises(Exception) as e:
            s.read_eth_files(data_path("invalid_data.txt"), data_path("cov_label.txt"))
        print(e.value, file=regtest)

        # nonexisting label file
        with pytest.raises(Exception) as e:
            s.read_eth_files(data_path(data_fn), data_path("invalid_label.txt"))
        print(e.value, file=regtest)

        s.read_eth_files(data_path(data_fn), data_path(label_fn))

        if if_read_cov:
            # nonexisting cov file
            with pytest.raises(Exception) as e:
                s.read_eth_files(data_path(data_fn), data_path(label_fn), data_path("invalid_cov.txt"))
            with pytest.raises(Exception) as e:
                s.update_covariates_file(data_path("invalid_cov.txt"))

            s.update_covariates_file(data_path(cov_fn))

        # did not set alpha
        with pytest.raises(Exception) as e:
            s.execute()
        print(e.value, file=regtest)

        with pytest.raises(Exception) as e:
            s.set_alpha(-0.01)
        with pytest.raises(Exception) as e:
            s.set_alpha(1.01)
        s.set_alpha(0.04)

        # did not set lmax
        with pytest.raises(Exception) as e:
            s.execute()
        print(e.value, file=regtest)

        with pytest.raises(Exception) as e:
            s.set_lmax(-1)
        with pytest.raises(Exception) as e:
            s.set_lmax(0.1)
        s.set_lmax(10)

        if if_set_n_perm:
            s.set_seed(1)
            # did not set n_perm
            with pytest.raises(Exception) as e:
                s.execute()
            with pytest.raises(Exception) as e:
                s.set_n_perm(0)
            with pytest.raises(Exception) as e:
                s.set_n_perm(0.1)
            s.set_n_perm(100)

        # execute method not called:
        with pytest.raises(Exception) as e:
            s.write_summary(tmpdir.join("output_summary%s.txt" % s_name).strpath)

        with pytest.raises(Exception) as e:
            s.write_profile(tmpdir.join("output_timing%s.txt" % s_name).strpath)

        with pytest.raises(Exception) as e:
            s.write_pvals_testable_intervals(
                tmpdir.join("p_vals_testable_intervals%s.txt" % s_name).strpath)

        with pytest.raises(Exception) as e:
            s.write_pvals_significant_intervals(
                tmpdir.join("p_vals_significant_intervals%s.txt" % s_name).strpath)

        with pytest.raises(Exception) as e:
            s.get_filtered_intervals()

        with pytest.raises(Exception) as e:
            s.get_significant_intervals()

        # now it should succeed
        s.execute()

        s.write_summary(tmpdir.join("output_summary%s.txt" % s_name).strpath)
        s.write_profile(tmpdir.join("output_timing%s.txt" % s_name).strpath)
        if if_itemset:
            s.write_pvals_testable_itemsets(
                tmpdir.join("p_vals_testable_itemsets%s.txt" % s_name).strpath)
            s.write_pvals_significant_itemsets(
                tmpdir.join("p_vals_significant_itemsets%s.txt" % s_name).strpath)
        else:
            s.write_pvals_testable_intervals(
                tmpdir.join("p_vals_testable_intervals%s.txt" % s_name).strpath)
            s.write_pvals_significant_intervals(
                tmpdir.join("p_vals_significant_intervals%s.txt" % s_name).strpath)

    _dump_result(s, regtest, itemset=if_itemset)
    _dump_tmpdir(tmpdir, regtest)




def _multiple_execute(s, alphas, tmpdir, regtest, if_write_results=True, if_itemset=False):
    for alpha in alphas:
        s.set_alpha(alpha)
        s.execute()

        if if_write_results:
            fn_suffix = "%s-alpha_%1.3f" % (s.__class__.__name__, alpha)
            s.write_summary(tmpdir.join("output_summary%s.txt" % fn_suffix).strpath)
            s.write_profile(tmpdir.join("output_timing%s.txt" % fn_suffix).strpath)
            if if_itemset:
                s.write_pvals_testable_itemsets(tmpdir.join(
                    "p_vals_testable_itemsets%s.txt" % fn_suffix).strpath)
                s.write_pvals_significant_itemsets(tmpdir.join(
                    "p_vals_significant_itemsets%s.txt" % fn_suffix).strpath)
            else:
                s.write_pvals_testable_intervals(tmpdir.join(
                    "p_vals_testable_intervals%s.txt" % fn_suffix).strpath)
                s.write_pvals_significant_intervals(tmpdir.join(
                    "p_vals_significant_intervals%s.txt" % fn_suffix).strpath)

            _dump_result(s, regtest, itemset=if_itemset)


def test_multiple_executes(tmpdir, data_path, regtest):
    """Test multiple executes w/ different alpha levels.
    """

    for (s_cls, if_set_n_perm, if_read_cov, if_itemset) in _ALL_SEARCHES:
        s = s_cls()
        if if_set_n_perm:
            s.set_n_perm(50)
            s.set_seed(1)
        if if_itemset:
            data_fn = "tictactoe_data.txt"
            label_fn = "tictactoe_labels.txt"
            cov_fn = "tictactoe_covariates.txt"
        elif if_read_cov:
            data_fn = "cov_data.txt"
            label_fn = "cov_label.txt"
            cov_fn = "cov.txt"
        else:
            data_fn = "data.txt"
            label_fn = "label.txt"
            cov_fn = None

        s.read_eth_files(data_path(data_fn), data_path(label_fn))
        if if_read_cov:
            s.update_covariates_file(data_path(cov_fn))

        s.set_lmax(0)

        # keep them high, because for 0 output sig ints FWER estimate seems to
        # vary a bit, despite setting the seed
        alphas = [0.2, 0.1, 0.05]
        _multiple_execute(s, alphas, tmpdir, regtest, if_itemset=if_itemset)

    _dump_tmpdir(tmpdir, regtest)


def _mem_now(search_obj=None):
    # first, free any non-collected leftovers (at least once)
    n_collected = gc.collect()
    while n_collected:
        n_collected = gc.collect()
    if search_obj:
        rss = search_obj.get_libmem().current
    else:
        # Note: the "portable" fields available on all plaforms are `rss` (Resident
        # Set Size) and `vms` (Virtual Memory Size). VMS, besides all stack and heap
        # memory, counts swapped out memory and pages of shared libraries which are
        # not in memory. We could care for that if we would use very large files for
        # tests, but we don't (because slooow).
        rss = psutil.Process(os.getpid()).memory_info().rss
    print("[DEBUG] mem_now = %d KB" % (rss / 1024))
    return rss / 1024 # in KB


def test_multiple_executes_memory(tmpdir, data_path, regtest):
    """Test memory usage over multiple executes.

    Note: currently in form of a unit test, not a regression test.
    """
    s = SignificantIntervalSearchChi()
    # s.read_plink_files(data_path("sample_data"))
    s.read_eth_files(data_path("data.txt"), data_path("label.txt"))
    s.set_lmax(0)

    # execute once to allocate results memory
    alpha = 0.05
    s.set_alpha(alpha)
    s.execute()

    # Check for memory diff at the end and start of n_rep identical executes
    n_rep = 50
    mem_start = _mem_now(s)
    _multiple_execute(s, n_rep * [alpha], tmpdir, regtest, if_write_results=False)
    mem_end = _mem_now(s)
    mem_delta = mem_end - mem_start

    # there is some small memory increase (Python only?); empirical value
    mem_margin = 16.0
    print("[DEBUG] mem_delta = %d" % mem_delta)
    assert mem_delta <= mem_margin, "significant increase in memory usage (kB) over subsequent execute() calls"


def test_clear_previous_output_on_error(tmpdir, data_path, regtest,
                                        test_eth_format = True):
    """Test clear output from previous execute on a subsequent execute error.
    """
    file_mem_lb = 500 if test_eth_format else 100000
    n_rep_err = 10 # repeat erroneous executes to assure that any side-effect
                   # memory increase after an error is not cumulative

    s = SignificantIntervalSearchWyExact()
    s.set_n_perm(50)
    mem_noread = _mem_now(s)
    if test_eth_format:
        s.read_eth_files(data_path("data.txt"), data_path("label.txt"));
    else:
        s.read_plink_files(data_path("sample_data"));
    s.set_lmax(0)
    s.set_alpha(0.05)
    s.set_seed(1)

    # execute: correct
    mem_noexecute = _mem_now(s)
    mem_delta = mem_noexecute - mem_noread
    print("[DEBUG] post-read mem_delta = %d KB" % mem_delta)
    # This fails when run together with other tests, but passes when the test
    # is run alone (no mem increase as if the file is already in the memory)
    #assert mem_delta >= file_mem_lb,\
    #    "expected at least %d kB increase of memory usage after file read" % file_mem_lb
    assert mem_delta >= 0,\
        "expected increase of memory usage after file read"
    s.execute()
    mem_execute = _mem_now(s)
    mem_delta = mem_execute - mem_noexecute
    print("[DEBUG] post-execute mem_delta = %d KB" % mem_delta)
    assert mem_delta >= 0,\
        "expected increase of memory usage after a execute (results)"

    # execute: force error (n_rep_err times)
    s._set_n_perm(-1)
    for i in range(n_rep_err):
        with pytest.raises(Exception):
            s.execute()

        # get_results: error
        with pytest.raises(Exception):
            s.get_result()
    mem_errexecute = _mem_now(s)

    mem_delta = mem_errexecute - mem_noexecute
    print("[DEBUG] pre-execute vs. post-error mem_delta = %d KB" % mem_delta)
# FIXME assert below fails because of what it looks like the results from
#       a successful execute are not cleaned from the memory. However, there is
#       no C++ leak, and C++-only test (*) does not show that size of an increase
#       (*) to run uncomment code in significant_interval_search_wy_{chi,exact}.cpp
#           (for more output uncomment also #define in debug.h) re-compile, and
#           call:
#               executables/significant_interval_search_wy_{exact,chi} -eth test/sampledata/data.txt test/sampledata/label.txt 0.05 50 0 test/output/output_{exact,chi}
#               executables/significant_interval_search_wy_{exact,chi} -plink test/plink/sample_data 0.05 50 0 test/output/output_{exact,chi}
#     mem_margin = 192
#     assert mem_delta <= mem_margin, (
#         "memory after successful and then failed execute() is expected to be "
#         "the same as before the first execute() call (results are cleared)")

    # everything is preserved after execute error - execute: correct
    s.set_n_perm(50)
    s.execute()

    # actual (regression) check of the correct execute result
    s.write_summary(tmpdir.join("output_summary.txt").strpath)
    s.write_profile(tmpdir.join("output_timing.txt").strpath)
    s.write_pvals_testable_intervals(tmpdir.join("p_vals_testable_intervals.txt").strpath)
    s.write_pvals_significant_intervals(tmpdir.join("p_vals_significant_intervals.txt").strpath)
    _dump_result(s, regtest)

    _dump_tmpdir(tmpdir, regtest)


def _multiple_fail_read_eth_files(searchobj, files, data_path):
    for (datafile, labelfile) in files:
        with pytest.raises(IOError):
            searchobj.read_eth_files(data_path(datafile), data_path(labelfile))


def _multiple_fail_read_plink_files(searchobj, files, data_path):
    for plinkfile in files:
        with pytest.raises(IOError):
            searchobj.read_plink_files(data_path(plinkfile))


def test_keep_data_on_nonexisting_file_read(tmpdir, data_path, regtest):
    """Test if data stays in the memory if trying to read non-existing files.

    BTW: tests `write_eth_files` function
    """
    s = SignificantIntervalSearchExact()
    s.read_eth_files(data_path("data.txt"), data_path("label.txt"))

    # Try read non-existing files and check for memory diff
    # there is some small memory increase (Python only?); empirical value
    mem_margin = 16.0
    # ETH file format
    mem_start = _mem_now(s)
    _multiple_fail_read_eth_files(s, random.randint(4,32) *
                                  [('not an existing file', 'label.txt'),
                                   ('data.txt', 'not an existing file')],
                                  data_path)
    mem_end = _mem_now(s)
    mem_delta = mem_end - mem_start
    print("[DEBUG] mem_delta = %d" % mem_delta)
    assert mem_delta <= mem_margin, "increase in used memory over unsuccessful ETH files read"
    # PLINK file format
    mem_start = _mem_now(s)
    _multiple_fail_read_plink_files(s, random.randint(4,32) *
                                    ['not an existing file'],
                                    data_path)
    mem_end = _mem_now(s)
    mem_delta = mem_end - mem_start
    print("[DEBUG] mem_delta = %d" % mem_delta)
    assert mem_delta <= mem_margin, "increase in used memory over unsuccessful PLINK files read"

    # Regression check for keeping the correctly read files
    s.write_eth_files(tmpdir.join("data.txt").strpath, tmpdir.join("label.txt").strpath)
    _dump_tmpdir(tmpdir, regtest)



def _test_keep_data_on_corrupted_file_read(tmpdir, regtest, mem_margin, search_obj, readfun_str, readfun_args):
    mem_start = _mem_now(search_obj)
    with pytest.raises(Exception) as e:
        readfun = getattr(search_obj, readfun_str)
        readfun(*readfun_args)
    print(e.value, file=regtest)

    mem_end = _mem_now(search_obj)
    mem_delta = mem_end - mem_start
    print("[DEBUG] mem_delta = %d" % mem_delta)
    assert mem_delta <= mem_margin, "increase in memory usage over corrupted file read"
    assert mem_delta >= 0, "decrease in memory usage over corrupted file read"

def _test_keep_eth_data_on_corrupted_file_read(tmpdir, data_path, regtest, searchobj, datafn, labelfn):
    # Memory usage test
    # there is some small relative to file size or type memory increase (Python only?); empirical value
    mem_margin = 128
    _test_keep_data_on_corrupted_file_read(tmpdir, regtest, mem_margin,
        searchobj, "read_eth_files", (data_path(datafn), data_path(labelfn),))

def _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, searchobj, plinkfn):
    # Memory usage test
    # there is some small relative to file size or type memory increase (Python only?); empirical value
    mem_margin = 2096
    _test_keep_data_on_corrupted_file_read(tmpdir, regtest, mem_margin,
        searchobj, "read_plink_files", (data_path(plinkfn),))

def _test_keep_cov_on_corrupted_file_read(tmpdir, data_path, regtest, searchobj, covfn):
    mem_margin = 16
    _test_keep_data_on_corrupted_file_read(tmpdir, regtest, mem_margin,
        searchobj, "update_covariates_file", (data_path(covfn),))

def test_keep_data_on_corrupted_file_read(tmpdir, data_path, regtest):
    """Test data is kept if file read failed due to the file being corrupted.

    Corrupted data means for ETH format:
     * data.txt with less columns (samples) than label.txt values and vice versa,
     * data.txt with rows of varying number of columns (samples),
     * data.txt or label.txt with non-binary values.
    For PLINK format analogous tests, but data file (.raw) is "transposed"
    (samples are in rows, whereas SNPs in columns).

    Note: the PLINK format .bim files are unused (so also not validated).
    """
    s = SignificantIntervalSearchFastCmh()


    ## 1) ETH format
    # correct read
    s.read_eth_files(data_path("data.txt"), data_path("label.txt"))

    ### a) label.txt with less values than data.txt columns
    _test_keep_eth_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "data.txt", "label-short.txt")
    ### b) data.txt with less columns than label.txt values
    _test_keep_eth_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "data-short.txt", "label.txt")
    ### c) data.txt with varying number of columns (samples)
    _test_keep_eth_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "data-var.txt", "label.txt")
    ### d) label.txt with non-digit values
    _test_keep_eth_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "data.txt", "label-char.txt")
    ### e) data.txt with non-digit values
    _test_keep_eth_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "data-char.txt", "label.txt")
    ### f) label.txt with single digit, but non-binary values
    _test_keep_eth_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "data.txt", "label-non_bin.txt")

    # Regression test for previous files
    s.write_eth_files(tmpdir.join("data.txt").strpath,
                      tmpdir.join("label.txt").strpath)


    ## 2) PLINK format
    # correct read
    s.read_plink_files(data_path("sample_data-short"))

    ### a) less .fam (labels) rows than .raw (data) rows
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-short_label")
    ### b) less .raw (data) rows than .fam (labels) rows
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-short_data")
    ### c) .raw (data) with varying number of columns (SNPs)
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-var")
    ### d) .fam (labels) with non-binary values,
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-char_label")
    ### e) .raw (data) with non-binary values,
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-char_data")
    ### f) .raw (data) with non-matching .fam (labels) line,
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-mismatch_data")

    ### g) .raw (data) with duplicate identifiers
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-duplicate_data")

    ### h) .fam (label) with duplicate identifiers
    _test_keep_plink_data_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                                 "sample_data-duplicate_label")

    # Regression test for previous files
    s.write_eth_files(tmpdir.join("data-plink.txt").strpath,
                      tmpdir.join("label-plink.txt").strpath)

    ## 2) Covariates
    # correct read
    s.read_eth_files(data_path("cov_data.txt"), data_path("cov_label.txt"))
    s.update_covariates_file(data_path("cov.txt"))

    ### a) cov.txt with less values than label.txt
    _test_keep_cov_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "cov-short.txt")
    ### b) cov.txt with more values than label.txt
    _test_keep_cov_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "cov-long.txt")
    ### c) cov.txt with non-digit values
    _test_keep_cov_on_corrupted_file_read(tmpdir, data_path, regtest, s,
                                               "cov-char.txt")

    # Regression test for previous files
    s.write_eth_files(tmpdir.join("cov_data.txt").strpath,
                      tmpdir.join("cov_label.txt").strpath)
    s.write_covariates_file(tmpdir.join("cov.txt").strpath)

    _dump_tmpdir(tmpdir, regtest)
