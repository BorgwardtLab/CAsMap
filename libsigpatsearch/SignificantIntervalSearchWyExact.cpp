/*
 * SignificantIntervalSearchWyExact.cpp
 *
 *  Created on: Aug 18, 2016
 *      Author: mabaker
 */

#include "SignificantIntervalSearchWyExact.h"

#include <stdlib.h>
#include <iostream>

#include "pval.h"



/* CONSTANT DEFINES */
#define NO_VERBOSE 1



using namespace std;

namespace SignificantPattern
{

SignificantIntervalSearchWyExact::SignificantIntervalSearchWyExact()
    : SignificantIntervalSearchWy()

{
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyExact()\n");
    #endif
    execute_constructor_wyexact();
}
SignificantIntervalSearchWyExact::~SignificantIntervalSearchWyExact()
{
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearchWyExact()\n");
    #endif
    execute_destructor_wyexact();
}



void SignificantIntervalSearchWyExact::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyExact::execute_constructor()\n");
    #endif
    execute_constructor_wyexact();
}
void SignificantIntervalSearchWyExact::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyExact::execute_destructor()\n");
    #endif
    execute_destructor_wyexact();
    super::execute_destructor();
}

void SignificantIntervalSearchWyExact::execute_constructor_wyexact() {
    loggamma_constructor();
    hypergeom_pvals_constructor();
}
void SignificantIntervalSearchWyExact::execute_destructor_wyexact(){
    hypergeom_pvals_destructor();
    loggamma_destructor();
}



void SignificantIntervalSearchWyExact::algorithm_init(){
    super::algorithm_init();

    // Initialise threshold value
    delta = ((double) n)/N; //$\psi(1)=\frac{n}{N}$
    loggamma_init();

    hypergeom_pvals_init();
}

void SignificantIntervalSearchWyExact::loggamma_init(){
    if (!loggamma) {
        loggamma = new double[N+1];
        loggamma_clear();
    }
}
void SignificantIntervalSearchWyExact::loggamma_clear(){
    longint x;
    // Initialise cache with appropriate values
    for(x=0;x<=N;x++) loggamma[x] = lgamma(x+1);//Gamma(x) = (x-1)!
    // Initialise log_inv_binom_N_n
    log_inv_binom_N_n = loggamma[n] + loggamma[N-n] - loggamma[N];
}
void SignificantIntervalSearchWyExact::loggamma_constructor(){
     loggamma = 0;
     log_inv_binom_N_n = 0;
}
void SignificantIntervalSearchWyExact::loggamma_destructor() {
    if (loggamma) delete [] loggamma;
    loggamma_constructor();
}

void SignificantIntervalSearchWyExact::hypergeom_pvals_init(){
    if (!hypergeom_pvals) {
        hypergeom_pvals = new double[n+1];
        hypergeom_pvals_clear();
    }
}
void SignificantIntervalSearchWyExact::hypergeom_pvals_clear(){
    std::fill_n(hypergeom_pvals, n+1, 0); // actual values computed in precompute_pvals
}
void SignificantIntervalSearchWyExact::hypergeom_pvals_constructor(){
    hypergeom_pvals = 0;
}
void SignificantIntervalSearchWyExact::hypergeom_pvals_destructor() {
    if (hypergeom_pvals) delete [] hypergeom_pvals;
    hypergeom_pvals_constructor();
}

void SignificantIntervalSearchWyExact::psi_clear(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWyExact::psi_clear()\n");
    #endif
    fisher_minpvals(N, n, N_over_2, psi);
}



/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

/* Evaluate Fisher's exact test on a table with margins x, n and N and cell count a. Note that n and N are defined as global variables.
 * The p-value is defined as a two-tailed p-value which adds up the probabilities of all tables less or equally likely to occur than the
 * one we observed
 */
double SignificantIntervalSearchWyExact::compute_pval(longint a, longint x){
    return fisher_pval(a, x, N, n, loggamma, log_inv_binom_N_n);
}

/* This function precomputes all Fisher exact test P-values for a contingency table with margins x,n,N that is,
 * all p-values p(a,x,n,N) for a in the range [max(0,n+x-N),min(x,n)]. The results will be stored in the array
 * hypergeom_pvals such that p(a,x,n,N)=hypergeom_pvals[a]. Note that values hypergeom_pvals[a] for a outside
 * [max(0,n+x-N),min(x,n)] are undefined and could contain garbage of previous hypotheses.
 * */
 void SignificantIntervalSearchWyExact::precompute_pvals(longint x){
    double pre_comp_xterms, pval, p_left, p_right;
    int a, a_min, a_max;

    // Compute the contribution of all terms depending on x but not on a
    pre_comp_xterms = loggamma[x] + loggamma[N-x];
    a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
    a_max = (x > n) ? n : x;//min(x,n)

    // Precompute the hypergeometric PDF in the range of interest
    for(a=a_min; a<=a_max; a++) hypergeom_pvals[a] = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a] + loggamma[n-a] + loggamma[x-a] + loggamma[(N-n)-(x-a)]));

    // The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
    // hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
    // that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. When a value is "accepted", we know its
    // respective p-value because due to the way we explore the hypergeometric pdf, there can be no other values larger than it. Therefore, everytime
    // a value is "accepted", we store the pvalue in hypergeom_pvals.
    // The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
    // that case is by "accepting" both values simultaneously.
    pval = 0;
    while(a_min<a_max){
        p_left = hypergeom_pvals[a_min]; p_right = hypergeom_pvals[a_max];
        if(p_left == p_right) { pval += (p_left+p_right); hypergeom_pvals[a_min++] = pval; hypergeom_pvals[a_max--] = pval; }
        else if(p_left < p_right){ pval += p_left; hypergeom_pvals[a_min++] = pval;}
        else{ pval += p_right; hypergeom_pvals[a_max--] = pval;}
    }
    // In this case a_min=a_max is the mode of the distribution and its p-value is 1 by definition
    if(a_min==a_max) hypergeom_pvals[a_max] = 1;
}

void SignificantIntervalSearchWyExact::process_first_layer_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    longint tau, queue_idx, i, j;
    unsigned char *X_tr_aux;
    unsigned char *Y_tr_perm_aux;
    double pval;
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
            // Precompute CDF and p-values of hypergeometric distribution with frequency x
            precompute_pvals(freq_par[tau]);
            // Compute cell-counts for all permutations
            for(i=0; i<N; i++){
                if(X_tr_aux[i]) for(j=0, Y_tr_perm_aux = Y_tr_perm[i]; j<n_perm; j++) a_cnt[j] += Y_tr_perm_aux[j];
            }
            // Check if we have a new minimum P-value for some of the permutations
            for(j=0; j<n_perm; j++){
                pval = hypergeom_pvals[a_cnt[j]]; a_cnt[j] = 0;
                if(pval < min_pval[j]){
                    // Increase FWER only if the previous minimum p-value is above current threshold
                    if((pval <= delta) && (min_pval[j] > delta)) FWER += ((double)1)/n_perm;
                    // Update minimum p-value
                    min_pval[j] = pval;
                }
            }
            update_threshold();
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

void SignificantIntervalSearchWyExact::process_intervals_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char **X_par = genotype_par.getMatrixPtr();
    longint tau, queue_idx, i, j;
    unsigned char *X_tr_aux, *X_par_aux;
    unsigned char *Y_tr_perm_aux;
    double pval;
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
            // Precompute CDF and p-values of hypergeometric distribution with frequency x
            precompute_pvals(freq_par[tau]);
            // Compute cell-counts for all permutations
            for(i=0; i<N; i++){
                if(X_par_aux[i]) for(j=0, Y_tr_perm_aux = Y_tr_perm[i]; j<n_perm; j++) a_cnt[j] += Y_tr_perm_aux[j];
            }
            // Check if we have a new minimum P-value for some of the permutations
            for(j=0; j<n_perm; j++){
                pval = hypergeom_pvals[a_cnt[j]]; a_cnt[j] = 0;
                if(pval < min_pval[j]){
                    // Increase FWER only if the previous minimum p-value is above current threshold
                    if((pval <= delta) && (min_pval[j] > delta)) FWER += ((double)1)/n_perm;
                    // Update minimum p-value
                    min_pval[j] = pval;
                }
            }
            update_threshold();
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
