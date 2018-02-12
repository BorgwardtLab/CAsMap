/*
 * SignificantIntervalSearchExact.h
 *
 *  Created on: Aug 18, 2016
 *      Author: mabaker
 */

#ifndef SIGNIFICANTINTERVALSEARCHWYEXACT_H_
#define SIGNIFICANTINTERVALSEARCHWYEXACT_H_

#include <string>
#include <iostream>
#include <fstream>

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */
#include <time.h> //Already included in original LCM source code above
#include <sys/time.h>
#include <sys/resource.h>
/* END OF LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */

#include "SignificantIntervalSearchWy.h"
#include "PValues.h"
#include "chi2.h"


namespace SignificantPattern
{

/// @todo
/// consider to multi-inherit from SignificantIntervalSearchExact and
/// SignificantIntervalSearchWy to factor out loggama_*, compute_pval,
/// psi_clear and so forth, but, again, beware of the deadly diamond of death
class SignificantIntervalSearchWyExact : public SignificantIntervalSearchWy
{

private:
    /**
     * super class pattern for code independence of changes in inheritance
     */
    typedef SignificantIntervalSearchWy super;

    double *loggamma;
    /**
     * Array for storing values of the PDF of the hypergeometric distribution
     * for fast p-value computation and eventually the p-values themselves (the
     * array is reused).
     */
    double *hypergeom_pvals;

protected:
    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    void execute_constructor_wyexact();
    void execute_destructor_wyexact();

    void loggamma_init();
    void loggamma_constructor();
    /**
     * Precompute values of log(x!) storing them in the array #loggamma
     */
    void loggamma_clear();
    void loggamma_destructor();

    /**
     * Allocate memory for #hypergeom_pvals (worst case memory requirement n+1).
     */
    void hypergeom_pvals_init();
    void hypergeom_pvals_constructor();
    /**
     * Pre-compute p-values of the hypergeometric distribution storing them in
     * #hypergeom_pvals.
     */
    void hypergeom_pvals_clear();
    void hypergeom_pvals_destructor();

    virtual void psi_clear() override;

    virtual void algorithm_init() override;
    // virtual void algorithm_end() override;


    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

    double compute_pval(longint a, longint x) override;

    void process_first_layer_threshold() override;
    void process_intervals_threshold() override;

    void precompute_pvals(longint x);

public:
    SignificantIntervalSearchWyExact();
    virtual ~SignificantIntervalSearchWyExact();

};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHWYEXACT_H_ */
