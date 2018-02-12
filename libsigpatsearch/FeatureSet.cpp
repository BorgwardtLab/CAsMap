/*
 * FeatureSet.cpp
 *
 *  Created on: Sep 8, 2016
 *      Author: mabaker
 */

#include <iostream> // hexfloat, defaultfloat
#include <sstream>

#include "FeatureSet.h"
#include "Exception.h"

namespace SignificantPattern
{

    FeatureSet::FeatureSet () : alphaVector(), pValueVector() {}
    FeatureSet::~FeatureSet () {}

    const std::string FeatureSet::COL_SEP = ",";
    const std::string FeatureSet::HEADER_PROPS = "p-value" + FeatureSet::COL_SEP + "a";

    void FeatureSet::addFeatureProps(longint alpha, double pValue)
    {
        alphaVector.push_back(alpha);
        pValueVector.push_back(pValue);
    }

    std::string const& FeatureSet::getHeaderProps() const {
        return HEADER_PROPS;
    }
    std::string const FeatureSet::getLineProps(size_t i) const {
        std::stringstream ss;
        ss << std::scientific << pValueVector[i] << COL_SEP
           << std::defaultfloat << alphaVector[i];
        return ss.str();
    }

    void FeatureSet::writeHeaderToFile(std::ofstream& file) const
    {
        file << getHeaderProps() << COL_SEP << getHeaderFeature() << std::endl;
    }
    void FeatureSet::writeLineToFile(std::ofstream& file, size_t i) const
    {
        file << getLineProps(i) << COL_SEP << getLineFeature(i) << std::endl;
    }

    void FeatureSet::writeToFile(const std::string& filename) const
    {
        std::ofstream file;
        file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
        try
        {
            file.open(filename.c_str());
        }
        catch (const std::ios_base::failure& e)
        {
            throw Exception("Failed opening " + filename + " for writing");
        }
        writeHeaderToFile(file);
        for (size_t i=0; i<getLength(); ++i) writeLineToFile(file, i);
        file.close();

    }





    ItemsetSet::ItemsetSet () : FeatureSet(), itemsetsVector() {}
    ItemsetSet::~ItemsetSet () {}

    const std::string ItemsetSet::ITEMS_SEP = " ";
    const std::string ItemsetSet::HEADER_FEATURE = "itemset";

    void ItemsetSet::addFeature(const std::vector<longint> itemset, longint alpha, double pValue)
    {
        addFeatureProps(alpha, pValue);
        itemsetsVector.push_back(itemset);
    }

    std::string const& ItemsetSet::getHeaderFeature() const {
        return HEADER_FEATURE;
    }
    std::string const ItemsetSet::getLineFeature(size_t i) const {
        std::stringstream ss;
        std::vector<longint> itemset = itemsetsVector[i];
        size_t n = itemset.size();
        for(size_t i = 0; i < n-1; i++)
            ss << itemset[i] << ITEMS_SEP;
        ss << itemset[n-1];
        return ss.str();
    }





    IntervalSet::IntervalSet () : FeatureSet(), startVector(), endVector() {}
    IntervalSet::~IntervalSet () {}

    const std::string IntervalSet::HEADER_FEATURE = "start" + IntervalSet::COL_SEP + "end";

    void IntervalSet::addFeature(longint start, longint end, longint alpha, double pValue)
    {
        addFeatureProps(alpha, pValue);
        startVector.push_back(start);
        endVector.push_back(end);
    }

    void IntervalSet::getLAndTauVectors(std::vector<longint>& lVector, std::vector<longint>& tauVector) const
    {
        for (size_t i=0; i<getLength(); ++i)
        {
            longint tau = startVector[i];
            longint l = endVector[i] - tau + 1;
            lVector.push_back(l);
            tauVector.push_back(tau);
        }
    }

    std::string const& IntervalSet::getHeaderFeature() const {
        return HEADER_FEATURE;
    }
    std::string const IntervalSet::getLineFeature(size_t i) const {
        std::stringstream ss;
        ss << startVector[i] << COL_SEP << endVector[i];
        return ss.str();
    }



    IntervalSetWithFreq::IntervalSetWithFreq() : IntervalSet(), xVector() {}
    IntervalSetWithFreq::~IntervalSetWithFreq() {}

    const std::string IntervalSetWithFreq::HEADER_PROPS_WITH_FREQ = IntervalSetWithFreq::HEADER_PROPS + IntervalSetWithFreq::COL_SEP + "x";

    void IntervalSetWithFreq::addFeature(longint start, longint end,
                                        longint alpha, double pValue) {
        addFeature(start, end, alpha, -1, pValue);
    }
    void IntervalSetWithFreq::addFeature(longint start, longint end,
                                        longint alpha, longint x,
                                        double pValue) {
        super::addFeature(start, end, alpha, pValue);
        xVector.push_back(x);
    }

    std::string const& IntervalSetWithFreq::getHeaderProps() const {
        return HEADER_PROPS_WITH_FREQ;
    }
    std::string const IntervalSetWithFreq::getLineProps(size_t i) const {
        std::stringstream ss;
        ss << super::getLineProps(i) << COL_SEP << xVector[i];
        return ss.str();
    }

} /* namespace SignificantPattern */
