from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "types.h":
    ctypedef long longint

cdef extern from "Summary.h" namespace "SignificantPattern":

    cdef cppclass Summary:
        longint getN()
        longint getn()
        longint getL()
        longint getm()
        longint getNumSignificantFeatures()
        longint getNumFeaturesProcessed()
        longint getNumFeatureSetsTotal()
        longint getL_max()
        double getDelta()
        double getDelta_opt()
        double getAlpha()
        void writeToFile(string filename) except +

    cdef cppclass SummaryInt(Summary):
        longint getMaxTestableIntervalLength()

    cdef cppclass SummaryFais(SummaryInt):
        longint getSl1()
        longint getSl2()
        longint getSu1()
        longint getSu2()

    cdef cppclass SummaryIset(Summary):
        pass

    cdef cppclass SummaryFacs(SummaryIset):
        longint getNumItemsetsClosedProcessed()

cdef extern from "Profiler.h" namespace "SignificantPattern":

    cdef cppclass Profiler:
        void writeToFile(string filename) except +
        size_t getPeakMemory()
        size_t getCurrMemory()


cdef extern from "FeatureSet.h" namespace "SignificantPattern":

    cdef cppclass IntervalSet:
        void writeToFile(string filename) except +

    cdef cppclass ItemsetSetWithOddsRatio:
        vector[vector[longint]] & getItemsetsVector()
        vector[double] & getScoreVector()
        vector[double] & getOddsRatioVector()
        vector[double] & getPValueVector()
        void writeToFile(string filename) except +


cdef extern from "FilterIntervals.h" namespace "SignificantPattern":

    cdef cppclass Interval:
        size_t getStart()
        size_t getEnd()
        double getScore()
        double getOddsRatio()
        double getPvalue()

    cdef cppclass FilterIntervals:
        vector[Interval] & getSigInts()
        void writeToFile(string filename) except +

    cdef cppclass SignificantIntervals:
        vector[Interval] & getSigInts()
        void writeToFile(string filename) except +

cdef extern from "SignificantFeaturesSearch.h" namespace "SignificantPattern":

    cdef cppclass SignificantFeaturesSearch:

        void readETHFiles(string, string, string) except +
        void writeETHFiles(string, string) except +
        void readPlinkFiles(string, string) except +
        void execute(double, longint) except +
        Profiler & getProfiler()

cdef extern from "SignificantIntervalSearch.h" namespace "SignificantPattern":

    cdef cppclass SignificantIntervalSearch(SignificantFeaturesSearch):

        FilterIntervals & getFilteredIntervals()
        SignificantIntervals & getSignificantIntervals()
        IntervalSet & getPValsTestableInts()
        IntervalSet & getPValsSigInts()


cdef extern from "SignificantIntervalSearchFais.h" namespace "SignificantPattern":

    cdef cppclass SignificantIntervalSearchFais(SignificantIntervalSearch):

        SummaryFais & getSummary()

cdef extern from "SignificantIntervalSearchExact.h" namespace "SignificantPattern":

    cdef cppclass SignificantIntervalSearchExact(SignificantIntervalSearchFais):

        SignificantIntervalSearchExact()


cdef extern from "SignificantIntervalSearchChi.h" namespace "SignificantPattern":

    cdef cppclass SignificantIntervalSearchChi(SignificantIntervalSearchFais):

        SignificantIntervalSearchChi()


cdef extern from "SignificantFeaturesSearchWithCovariates.h" namespace "SignificantPattern":

    cdef cppclass SignificantFeaturesSearchWithCovariates(SignificantFeaturesSearch):

        void readETHFilesWithCovariates(string, string, string, bool, string) except +
        void writeETHFilesWithCovariates(string, string, string) except +
        void readPlinkFilesWithCovariates(string, string, bool, string) except +
        void readCovariatesFile(string) except +
        void writeCovariatesFile(string) except +

cdef extern from "SignificantIntervalSearchFastCmh.h" namespace "SignificantPattern":

    cdef cppclass SignificantIntervalSearchFastCmh(SignificantIntervalSearch,SignificantFeaturesSearchWithCovariates):

        SummaryInt & getSummary()

        SignificantIntervalSearchFastCmh()

cdef extern from "SignificantItemsetSearch.h" namespace "SignificantPattern":

    cdef cppclass SignificantItemsetSearch:

        ItemsetSetWithOddsRatio & getPValsTestableIsets()
        ItemsetSetWithOddsRatio & getPValsSigIsets()

cdef extern from "SignificantItemsetSearchFacs.h" namespace "SignificantPattern":

    cdef cppclass SignificantItemsetSearchFacs(SignificantItemsetSearch,SignificantFeaturesSearchWithCovariates):

        SummaryFacs & getSummary()

        SignificantItemsetSearchFacs()
