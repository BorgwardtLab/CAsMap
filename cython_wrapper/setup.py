# encoding: utf-8
from __future__ import print_function, division, absolute_import
import sys

from setuptools import setup, Extension

DEBUG = False
if "--debug" in sys.argv:
    DEBUG = True
    sys.argv.remove("--debug")

# meta data of package #######################################

license = "GPL (>= 2)"

author = "Matthew Baker, Dean Bodenham, Felipe Llinares Lopez, Mikolaj Rybinski, Uwe Schmitt"
author_email = "mikolaj.rybinski@id.ethz.ch"
url = "https://sissource.ethz.ch/sispub/significant_interval_search_rewrite"

description = """Significant Pattern Search for Genome Analysis"""

long_description = """
A method based on significant pattern mining that is used to detect intervals in
binary genotype data which are significantly associated with a particular
phenotype.
"""

version = "0.3"


# internal stuff #############################################
package = "sigpatsearch"

requires = [
    'cython',
]

tests_require = [
    'psutil',
    'pytest',
    'pytest-regtest',
]


import glob
sources = [f.strip() for f in glob.glob("sigpatsearch/*.cpp")]

extra_compile_args = ['-std=c++11',] if not DEBUG \
    else ['-O0', '-g3', '-DDEBUG', '-Wall', '-pedantic', '-std=c++11',]

module = Extension(
    'sigpatsearch.wrapper',
    sources = sources,
    language = "c++",
    extra_compile_args = extra_compile_args,
    extra_link_args = extra_compile_args,
)

setup(
    name="sigpatsearch",
    version=version,
    license=license,
    description=description,
    long_description=long_description,
    url=url,
    author=author,
    author_email=author_email,
    packages=[package],
    include_package_data=True,
    install_requires=requires,
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=tests_require,
    test_suite="pytest",
    ext_modules=[module],
)
