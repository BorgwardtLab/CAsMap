/*
 * SignificantIntervalSearchWyChi.cpp
 *
 *  Created on: Aug 18, 2016
 *      Author: mabaker
 */

#include "SignificantIntervalSearchWyChi.h"

#include <stdlib.h>
#include <iostream>

#include "pval.h"



/* CONSTANT DEFINES */
#define NO_VERBOSE 1



using namespace std;

namespace SignificantPattern
{


SignificantIntervalSearchWyChi::SignificantIntervalSearchWyChi()
    : SignificantIntervalSearchWy()
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyChi()\n");
    #endif
    execute_constructor_wychi();
}
SignificantIntervalSearchWyChi::~SignificantIntervalSearchWyChi() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearchWyChi()\n");
    #endif
    execute_destructor_wychi();
}


void SignificantIntervalSearchWyChi::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyChi::execute_constructor()\n");
    #endif
    execute_constructor_wychi();
}
void SignificantIntervalSearchWyChi::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyChi::execute_destructor()\n");
    #endif
    execute_destructor_wychi();
    super::execute_destructor();
}


void SignificantIntervalSearchWyChi::execute_constructor_wychi() {
    class_ratio = 0; class_ratio_bin = 0;
}
void SignificantIntervalSearchWyChi::execute_destructor_wychi(){
}


void SignificantIntervalSearchWyChi::algorithm_init(){
    // Precompute constants for psi_init (called from algorithm_init)
    class_ratio = ((double)n)/N; class_ratio_bin = class_ratio*(1-class_ratio);
    super::algorithm_init();
    // Initialise threshold value
    delta = psi[1];
}
void SignificantIntervalSearchWyChi::algorithm_end(){
    super::algorithm_end();
}


void SignificantIntervalSearchWyChi::psi_clear(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyChi::psi_clear()\n");
    #endif
    chi2_minpvals(N, n, N_over_2, class_ratio, class_ratio_bin, psi);
}


/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
double SignificantIntervalSearchWyChi::compute_pval(longint a, longint x){
    return chi2_pval(a, x, N, n, class_ratio_bin);
}



void SignificantIntervalSearchWyChi::process_first_layer_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    longint tau, queue_idx, i, j;
    unsigned char *X_tr_aux;
    unsigned char *Y_tr_perm_aux;
    double Tval, pval, aux, num_precomp, den_precomp;
    // Process each length 1 interval
    for(tau=0; tau<L; tau++){
        n_featuresets_processed++;
        // Compute number of 1s in the interval
        X_tr_aux = X_tr[tau];
        for(j=0; j<N; j++) freq_par[tau] += X_tr_aux[j];
        // If the interval is testable...
        // Update frequency-buckets and number of testable intervals
        #ifndef NO_SINGLE_FEATURES
        if(istestable_int(tau)){
            // Precompute common parts of Chi2 statistic
            aux = ((double)freq_par[tau])/N; num_precomp = -n*aux; den_precomp = freq_par[tau]*(1-aux)*class_ratio_bin;
            if(den_precomp>0){//If den_precomp==0, then pval=1 and the min_pvals will never change
                // Compute cell-counts for all permutations
                for(i=0; i<N; i++){
                    if(X_tr_aux[i]) for(j=0, Y_tr_perm_aux = Y_tr_perm[i]; j<n_perm; j++) a_cnt[j] += Y_tr_perm_aux[j];
                }
                // Check if we have a new minimum P-value for some of the permutations
                for(j=0; j<n_perm; j++){
                    Tval = a_cnt[j]+num_precomp; Tval *= Tval; Tval /= den_precomp;
                    pval = Chi2_sf(Tval,1); a_cnt[j] = 0;
                    if(pval < min_pval[j]){
                        // Increase FWER only if the previous minimum p-value is above current threshold
                        if((pval <= delta) && (min_pval[j] > delta)) FWER += ((double)1)/n_perm;
                        // Update minimum p-value
                        min_pval[j] = pval;
                    }
                }
                update_threshold();
            }
        }
        #endif
        // If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
        // then do NOT append the left-child to the testable queue (i.e. prune it)
        if((tau==0) || isprunable_int(tau) || isprunable_int(tau-1)) continue;
        // Compute index of array position to be used, wrapping around if necessary
        queue_idx = testable_queue_front + testable_queue_length;
        queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
        // Actually append children to testable queue
        testable_queue[queue_idx] = tau-1;
        // Update queue length
        testable_queue_length++;
    }
}

void SignificantIntervalSearchWyChi::process_intervals_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char **X_par = genotype_par.getMatrixPtr();
    longint tau, queue_idx, i, j;
    unsigned char *X_tr_aux, *X_par_aux;
    unsigned char *Y_tr_perm_aux;
    double Tval, pval, aux, num_precomp, den_precomp;
    // While testable-interval queue is not empty, continue to process intervals
    while(testable_queue_length){
        // Pop a testable interval from the queue
        tau = testable_queue[testable_queue_front];
        testable_queue_front = (testable_queue_front<(L-1)) ? testable_queue_front + 1: 0;
        testable_queue_length--;
        // Check if we have started processing a new layer by detecting non-monotonicity in tau
        if(tau < last_tau) {
            l++;
            #ifndef NO_VERBOSE
            printf("\tProcessing layer %lld...\n",l+1);
            #endif
        }
        if((L_max>0) && ((l+1) > L_max)) {
            #ifndef NO_VERBOSE
            printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n",l+1);
            #endif
            break;
        }
        last_tau = tau;
        // Check any of the two parents is prunable, stop processing. Notice that this check is necessary
        // even if the current interval was appended to the testable queue, because the threshold and
        // testability regions might have been modified between the time in which the current interval
        // was appended to the queue and the time in which it is being processed
        if(isprunable_int(tau) || isprunable_int(tau+1)) continue;
        n_featuresets_processed++;
        // Compute OR and frequency of the interval
        X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
        for(j=0; j<N; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; freq_par[tau]++;}
        // If the interval is testable, increase counter of testable items and frequency-buckets and
        // check if the corrected significance threshold must be reduced
        if(istestable_int(tau)){
            // Precompute common parts of Chi2 statistic
            aux = ((double)freq_par[tau])/N; num_precomp = -n*aux; den_precomp = freq_par[tau]*(1-aux)*class_ratio_bin;
            if(den_precomp>0){//If den_precomp==0, then pval=1 and the min_pvals will never change
                // Compute cell-counts for all permutations
                for(i=0; i<N; i++){
                    if(X_par_aux[i]) for(j=0, Y_tr_perm_aux = Y_tr_perm[i]; j<n_perm; j++) a_cnt[j] += Y_tr_perm_aux[j];
                }
                // Check if we have a new minimum P-value for some of the permutations
                for(j=0; j<n_perm; j++){
                    Tval = a_cnt[j]+num_precomp; Tval *= Tval; Tval /= den_precomp;
                    pval = Chi2_sf(Tval,1); a_cnt[j] = 0;
                    if(pval < min_pval[j]){
                        // Increase FWER only if the previous minimum p-value is above current threshold
                        if((pval <= delta) && (min_pval[j] > delta)) FWER += ((double)1)/n_perm;
                        // Update minimum p-value
                        min_pval[j] = pval;
                    }
                }
                update_threshold();
            }
        }
        // If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
        // then do NOT append the left-child to the testable queue (i.e. prune it)
        if((tau==0) || isprunable_int(tau) || isprunable_int(tau-1)) continue;
        // Compute index of array position to be used, wrapping around if necessary
        queue_idx = testable_queue_front + testable_queue_length;
        queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
        // Actually append children to testable queue
        testable_queue[queue_idx] = tau-1;
        // Update queue length
        testable_queue_length++;
    }
}



} /* namespace SignificantPattern */
