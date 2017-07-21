/*
 * SignificantIntervalSearchExact.h
 *
 *  Created on: Aug 18, 2016
 *      Author: mabaker
 */

#ifndef SIGNIFICANTINTERVALSEARCHWYCHI_H_
#define SIGNIFICANTINTERVALSEARCHWYCHI_H_

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
#include "Summary.h"

namespace SignificantPattern
{

/// @todo
/// consider to multi-inherit from SignificantIntervalSearchChi and
/// SignificantIntervalSearchWy to factor out compute_pval, psi_clear and so
/// forth, but, again, beware of the deadly diamond of death
class SignificantIntervalSearchWyChi : public SignificantIntervalSearchWy
{
private:
    // super class pattern for code independence of changes in inheritance
    typedef SignificantIntervalSearchWy super;


    double class_ratio, class_ratio_bin;

    /* -------------------------------- INITIALISATION AND TERMINATION METHODS ----------------------------------------- */
    virtual void execute_constructor() override;
    virtual void execute_destructor() override;

    void execute_constructor_wychi();
    void execute_destructor_wychi();

    virtual void algorithm_init() override;
    virtual void algorithm_end() override;

    virtual void psi_clear() override;


    /* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

    double compute_pval(longint a, longint x) override;

    void process_first_layer_threshold() override;
    void process_intervals_threshold() override;

public:
    SignificantIntervalSearchWyChi();
    virtual ~SignificantIntervalSearchWyChi();

};




} /* namespace SignificantPattern */

#endif /* SIGNIFICANTINTERVALSEARCHWYCHI_H_ */
