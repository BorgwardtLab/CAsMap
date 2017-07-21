/*
 * SignificantIntervalSearchFastCmh.h
 *
 *  Created on: 2 Mar 2017
 *      Author: mikolajr
 */

#ifndef SIGNIFICANTINTERVALSEARCHFASTCMH_H_
#define SIGNIFICANTINTERVALSEARCHFASTCMH_H_

#include "SignificantIntervalSearch.h"
#include "SignificantFeaturesSearchTaroneCmh.h"

namespace SignificantPattern
{

class SignificantIntervalSearchFastCmh : public SignificantIntervalSearch,
                                         public SignificantFeaturesSearchTaroneCmh
{
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantFeaturesSearchTaroneCmh super_cov;
    typedef SignificantIntervalSearch super_feat;

    IntervalSet pValsTestableInts;
    IntervalSet pValsSigInts;


    void execute_constructor_fastcmh();
    void execute_destructor_fastcmh();

protected:
    //OPT: Record frequencies for each class freq_par_cov[tau] (define and use
    //     IntervalSetWithFreqInClasses)
    inline void saveSignificantInterval(double pval, longint tau, longint l, longint a) override {
        pValsSigInts.addFeature(tau, tau+l, a, pval);
    }
    inline void saveTestableInterval(double pval, longint tau, longint l, longint a) override {
        pValsTestableInts.addFeature(tau, tau+l, a, pval);
    }

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    inline virtual void algorithm_init() override {
        super_cov::algorithm_init();
        super_feat::algorithm_init();
    };
    inline virtual void algorithm_end() override {
        super_feat::algorithm_end();
        super_cov::algorithm_end();
    }

    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
    inline double compute_interval_pval(longint a, longint tau) override {
        return compute_pval(a,freq_par_cov[tau]);
    };

    inline bool istestable_int(longint tau) override {
        return istestable_freqcov(freq_par_cov[tau]);
    }
    inline bool isprunable_int(longint tau) override {
        return isprunable_freqcov(freq_par_cov[tau]);
    };

    void process_first_layer_pvalues() override;
    void process_intervals_pvalues() override;

    void process_first_layer_threshold() override;
    void process_intervals_threshold() override;

public:
    SignificantIntervalSearchFastCmh();
    virtual ~SignificantIntervalSearchFastCmh();

    inline IntervalSet const& getPValsTestableInts() const override {
        return pValsTestableInts;
    }
    inline IntervalSet const& getPValsSigInts() const override {
        return pValsSigInts;
    }
};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHFASTCMH_H_ */
