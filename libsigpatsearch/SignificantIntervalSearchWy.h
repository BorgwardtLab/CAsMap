/*
 * SignificantIntervalSearchWy.h
 *
 *  Created on: Dec 20, 2016
 *      Author: mikolajr
 */

#ifndef LIBSIGINTERVALSEARCH_SIGNIFICANTINTERVALSEARCHWY_H_
#define LIBSIGINTERVALSEARCH_SIGNIFICANTINTERVALSEARCHWY_H_

#include "SignificantIntervalSearchFais.h"
#include "Summary.h"



namespace SignificantPattern
{

class SignificantIntervalSearchWy: public SignificantIntervalSearchFais
{

private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantIntervalSearchFais super;

    SummaryWy summary;

protected:
    inline virtual SummaryWy& getSummaryRef() override { return summary; }
    inline virtual void summary_constructor() override { summary = SummaryWy(); }

    // Variables for Westfall-Young permutation
    // Number of permutations
    int n_perm;

    // Matrix of size (# permutations) x (# non-empty transactions)
    unsigned char **Y_tr_perm;

    int* perm_array;
    std::string permutationsFilename;
    //Random generator seed for reproducibility
    int seed;


    // Current FWER
    double FWER;
    // FWER at corrected significance threshold
    double FWER_opt;
    // Minimum P-value for each permutation
    double *min_pval;

    // Cell-count counter (table entry (X=1,Y=1))
    // While finding significance threshold J-dimensional vector (one cell count per permutation)
    longint *a_cnt;

    /* --------------- INITIALISATION AND TERMINATION METHODS --------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    void execute_constructor_wy();
    void execute_destructor_wy();

    virtual void algorithm_init() override;
    virtual void algorithm_end() override;

    void perm_destructor();
    void perm_init();

    int rand_int(int x);
    void randperm(unsigned char *buffer, unsigned char *src, int n);
    void perm_array_init(int n);
    void perm_array_destructor();

    /* ------- FUNCTIONS TO FIND THE CORRECTED SIGNIFICANCE THRESHOLD ------- */
    void compute_corrected_significance_threshold() override;

    void decrease_threshold() override;
    inline void update_threshold() override {
        while(FWER > alpha) decrease_threshold();
    }


    // These have to be overwritten for WY methods
    virtual void process_first_layer_threshold() override = 0;
    virtual void process_intervals_threshold() override = 0;

public:
    SignificantIntervalSearchWy();
    virtual ~SignificantIntervalSearchWy();

    inline longint getNPerm() { return n_perm;}

    inline void setNPerm(longint val) { if (n_perm != val) perm_destructor(); n_perm = val; }
    inline void setPermutationsFilename(const std::string& filename) { permutationsFilename = filename;}
    inline virtual SummaryWy const& getSummary() const override { return summary; }
    void setSeed(unsigned s) { seed = s; };

};



}  /* namespace SignificantPattern */

#endif /* LIBSIGINTERVALSEARCH_SIGNIFICANTINTERVALSEARCHWY_H_ */
