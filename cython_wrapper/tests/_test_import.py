# encoding: utf-8
from __future__ import print_function, division, absolute_import

import glob
import hashlib
import os
from os.path import join, abspath, dirname

import pytest


@pytest.fixture
def path():
    here = dirname(abspath(__file__))

    def _(file_name):
        return join(here, file_name)
    return _


def _dump_outputs(tmpdir, regtest):

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


def test_import(path, tmpdir, regtest):
    from sigpatsearch import SignificantIntervalSearchExact
    search = SignificantIntervalSearchExact()
    search.read_eth_files(path("data.txt"), path("label.txt"))
    search.write_eth_files(path("data_copy.txt"), path("label_copy.txt"))
    search.set_alpha(0.4)
    search.set_lmax(0)
    search.execute()
    search.write_summary(tmpdir.join("summary.txt").strpath)
    search.write_profile(tmpdir.join("profile.txt").strpath)
    search.write_filtered_intervals(tmpdir.join("filtered_intervals.txt").strpath)
    search.write_pvals_testable_intervals(tmpdir.join("pvals_testable_intervals.txt").strpath)
    search.write_pvals_significant_intervals(tmpdir.join("pvals_significant_intervals.txt").strpath)
    search.get_significant_intervals()
    search.get_filtered_intervals()
    search.get_summary()
    search.get_result()
    search.get_libmem()
    _dump_outputs(tmpdir, regtest)
