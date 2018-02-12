/*
 * FilterIntervals.h
 *
 *  Created on: Sep 12, 2016
 *      Author: mabaker
 */

#ifndef FILTERINTERVALS_H_
#define FILTERINTERVALS_H_

//filterIntervals.h
//      read csv file containing all significant intervals, group these into
//      overlapping clusters, and return the most significant interval per cluster
//
//Dean Bodenham June 2016

#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<iomanip>
#include<numeric>

#include "types.h"

using std::vector;
using std::size_t;
using std::string;
using std::stringstream;

const double DEFAULT_PVALUE = 1.0;
const size_t DEFAULT_START = 0;
const size_t DEFAULT_END = 0;

namespace SignificantPattern
{
    //--------------------------------------------------------------------------//
    //--------------------------------------------------------------------------//
    class Interval{
        public:
            size_t getStart() const;
            size_t getEnd() const;
            double getPvalue() const;
            size_t getLength() const;
            void setStart(size_t);
            void setEnd(size_t);
            void setEnd(size_t, size_t);
            void setPvalue(double);
            bool overlaps(size_t, size_t) const;
            void printInterval() const;
            static size_t computeEnd(size_t tau_, size_t l_);
        private:
            size_t start;
            size_t end;
            double pvalue;
    };
    //--------------------------------------------------------------------------//

    class FilterIntervals
    {
    private:
        vector<Interval> sigInts;
        vector<Interval> getMinPvalueIntervalPerCluster(vector<size_t>& tau, vector<size_t>& l, vector<double>& pvalue, const vector<int>& label);
        vector<int> getClusterLabelsForIntervals(const vector<size_t>& tau, const vector<size_t>& l, const vector<Interval>& cluster);
        vector<Interval> getClusters(vector<size_t>& v_tau, vector<size_t>& v_l);
        vector<bool> getClusterIndicatorVector(vector<size_t>& v_tau, vector<size_t>& v_l);
        void makeIntervalTrue(vector<bool>& v, const size_t tau, const size_t l);
    public:
        FilterIntervals ();
        FilterIntervals (const FilterIntervals& other);
        virtual FilterIntervals& operator=(const FilterIntervals& other);
        virtual ~FilterIntervals ();
        void cpp_filterIntervalsFromMemory(vector<longint> ll_tau,
                                                         vector<longint> ll_l,
                                                    vector<double> pvalue);
        void writeToFile(const std::string& filename);
        inline vector<Interval>& getSigInts() { return sigInts; }
    };

    class SignificantIntervals
    {
    private:
        vector<Interval> sigInts;

    public:
        SignificantIntervals();
        SignificantIntervals (const SignificantIntervals& other);
        virtual SignificantIntervals& operator=(const SignificantIntervals& other);
        virtual ~SignificantIntervals();
        void cpp_intervalsFromMemory(vector<longint> ll_tau,
                                     vector<longint> ll_l,
                                     vector<double> pvalue);
        void writeToFile(const std::string& filename);
        inline vector<Interval>& getSigInts() { return sigInts; }
    };

} /* namespace SignificantPattern */

#endif /* FILTERINTERVALS_H_ */
