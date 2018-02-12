/*
 * SignificantIntervalSearchWy.cpp
 *
 *  Created on: Dec 20, 2016
 *      Author: trybik
 */

#include "SignificantIntervalSearchWy.h"

#include <math.h> //floor
#include <stdlib.h>
#include <vector>
#include<algorithm> //sort

#include "Exception.h"
//#include "double_comp.h" //doublecomp for C99 qsort



/* CONSTANT DEFINES */
#define NO_VERBOSE 1



namespace SignificantPattern
{



SignificantIntervalSearchWy::SignificantIntervalSearchWy()
    : SignificantIntervalSearchFais(),
      n_perm(0), permutationsFilename(""), seed(-1)
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy()\n");
    #endif
    execute_constructor_wy();
}
SignificantIntervalSearchWy::~SignificantIntervalSearchWy() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearchWy()\n");
    #endif
    execute_destructor_wy();
}



void SignificantIntervalSearchWy::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy::execute_constructor()\n");
    #endif
    execute_constructor_wy();
}
void SignificantIntervalSearchWy::execute_destructor() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy::execute_destructor()\n");
    #endif
    execute_destructor_wy();
    super::execute_destructor();
}


void SignificantIntervalSearchWy::execute_constructor_wy() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy::execute_constructor_wy()\n");
    #endif
    // TODO keep Y_tr_perm or perm_array mem in between calls if, respectively,
    //      n_perm and N or N have not changed
    Y_tr_perm = 0; perm_array = 0;

    FWER = 0; FWER_opt = 0;
    min_pval = 0; a_cnt = 0;
}
void SignificantIntervalSearchWy::execute_destructor_wy() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy::execute_destructor_wy()\n");
    #endif
    if (min_pval) 
        delete [] min_pval; 
    min_pval = 0;
    if (a_cnt) 
        delete [] a_cnt; 
    a_cnt = 0;

    perm_destructor();
}


void SignificantIntervalSearchWy::algorithm_init(){
    super::algorithm_init();

    unsigned char *Y_tr = phenotype.getVectorPtr();

    // Initialise Westfall-Young permutation-related variables
    // Save core constants as global variables
    perm_init();

    /* Compute the random permutations of the class labels */
    // Initialize random number generator
    if (seed == -1) seed = time(NULL);
    srand(seed);
    // Do the n_perm permutations themselves, storing them in labels_perm and
    // saving them to the file f_perm if necessary
    vector<unsigned char> perm_buf(N);
    unsigned char *perm_buf_start = &perm_buf[0];
    for(int j=0;j<n_perm;j++) {
        //TODO longint N , whereas randperm( , , int)
        randperm(perm_buf_start, Y_tr, N);
        // Dump contents of buffer into destination, skipping values corresponding to empty observations/transactions
        for(int i=0;i<N;i++) Y_tr_perm[i][j] = perm_buf_start[i];
    }

    // Allocate memory for minimum p-values, raising error if it fails
    min_pval = new double[n_perm];
    // Initialise all p-values to 1
    std::fill_n(min_pval, n_perm, 1);

    // Allocate memory for cell counts, raising an error if it fails
    a_cnt = new longint[n_perm];
    std::fill_n(a_cnt, n_perm, 0);
}
void SignificantIntervalSearchWy::algorithm_end(){
    super::algorithm_end();

    SummaryWy& summary = getSummaryRef();
    summary.setJ(n_perm);

    summary.setFWER(FWER);
    summary.setFWER_opt(FWER_opt);
    summary.setMin_pval(min_pval);
}



void SignificantIntervalSearchWy::perm_init() {
    if (n_perm <= 0)
        throw Exception("Invalid parameter value: set a positive number of permutations via setNPerm() call");

    perm_destructor();

    // If the generated permutation is to be stored in a file, initialise memory for perm_array
    if(permutationsFilename != ""){
        perm_array_init(N);
    }

    /* Allocate memory for the matrix of permuted labels */
    // First allocate memory for the n_perm row pointers
    Y_tr_perm = new unsigned char *[N];

    // Now allocate memory for a contiguous block of n_perm*(# non-empty transactions) chars
    Y_tr_perm[0] = new unsigned char[n_perm*N];
    // And make each row pointer point to the appropriate position (Y_tr_perm[0]
    // is already correctly set)
    for(longint i=1;i<N;i++) Y_tr_perm[i] = Y_tr_perm[0] + i*n_perm;
}
void SignificantIntervalSearchWy::perm_destructor() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy::perm_destructor()\n\tY_tr_perm = %p\n", (void *) Y_tr_perm);
    #endif
    if (Y_tr_perm)
    {
        #ifdef DEBUG
        fprintf(stderr, "\tY_tr_perm[0] = %p\n", (void *) Y_tr_perm[0]);
        #endif
        if (Y_tr_perm[0]) delete [] Y_tr_perm[0];
        delete [] Y_tr_perm; Y_tr_perm = 0;
    }
    perm_array_destructor();
}



/* -----------------FUNCTIONS TO SAMPLE RANDOM INTEGERS AND GENERATE RANDOM PERMUTATIONS--------------------------- */
void SignificantIntervalSearchWy::perm_array_init(int n) {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy::perm_array_init()\n\tperm_array = %p\n", (void *) perm_array);
    #endif
    //TODO cleanup & allocate only if N has changed (on a different file read)
    perm_array_destructor();
    perm_array = new int[n];
    // Initialise indices of the permutation to the identity permutation
    for (int i =0; i < n; i++) perm_array[i] = i; // C++11: std::iota(perm_array, perm_array+n, 0);
}
void SignificantIntervalSearchWy::perm_array_destructor() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchWy::perm_array_destructor()\n\tperm_array = %p\n", (void *) perm_array);
    #endif
    if (perm_array) 
        delete [] perm_array; 
    perm_array = 0;
}

int SignificantIntervalSearchWy::rand_int(int x){
/**
 * Sample a random integer uniformly distributed the range [0,x)
 *
 * Default random number generator samples integers uniformly in the range [0,RAND_MAX)
 * To sample in the range [0,x) with x < RAND_MAX, one can think of sampling an integer
 * rnd from [0,RAND_MAX) and then returning rnd % x. However, this may generate a non-uniform
 * distribution unless RAND_MAX is a multiple of x. To fix that, we combine that basic idea
 * with rejection sampling, rejecting any value rnd greater or equal than RAND_MAX - RAND_MAX *x.
 * This is the same as ensuring that the "effective" RAND_MAX is an exact multiple of x, so that
 * returning rnd % x leads to a uniform sampling scheme in the range [0,x)
 * The seed of the random number generator should have been initialised externally
 *
 */
    int rnd;
    int limit = RAND_MAX - RAND_MAX % x;

    do{
        rnd = rand();
    }while(rnd >= limit);
    return rnd % x;
}

void SignificantIntervalSearchWy::randperm(unsigned char *buffer, unsigned char *src, int n){
    int i,j; // Variables for looping and swapping
    char tmp; // Temp int for swapping
    //Array to store the permutation and temp int for swapping (only needed if permutationsFilename is not empty)
    int tmp_int;

    // First of all, copy the original array in the buffer
    std::copy(src, src+n, buffer);

    // Fisher-Yates algorithm
    for(i = n-1; i > 0; i--){
        // Sample a random integer in [0,i]
        j = rand_int(i + 1);
        // Swap dest[j] and dest[i]
        tmp = buffer[j];
        buffer[j] = buffer[i];
        buffer[i] = tmp;
        // If the generated permutation has to be saved, keep
        // track of the indices as well
        if (perm_array) {
            tmp_int = perm_array[j];
            perm_array[j] = perm_array[i];
            perm_array[i] = tmp_int;
        }
    }
    // If the generated permutation is to be stored in a file, write perm_array to the stream
    /*if(f_perm){
        // Write the first l-1 indices with comma as a delimiter
        for(i=0;i<(l-1);i++) fprintf(f_perm,"%d,",perm_array[i]);
        // For the last index, change the comma by a newline char
        fprintf(f_perm,"%d\n",perm_array[i]);
        // Free allocated memory
        free(perm_array);
    }*/
}

/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
/* Wrapper function that encapsulates the functionality required to find the corrected significance threshold */
void SignificantIntervalSearchWy::compute_corrected_significance_threshold(){
    FWER = 0; // Important: before super call
    // compute non-WY threshold
    super::compute_corrected_significance_threshold();
    // Set final corrected significance threshold
    // Sort p-values
    std::sort(min_pval, min_pval+n_perm); //C99: qsort(min_pval,n_perm,sizeof(double),doublecomp);
    // Tentative index to corrected significance threshold
    longint idx_max = (int)floor(alpha*n_perm)-1;
    // override non-WY threshold
    delta_opt = min_pval[idx_max];
    // Check and correct (if necessary) boundary cases
    if(delta_opt==min_pval[idx_max+1]){
        while(min_pval[--idx_max]==delta_opt);
        delta_opt = min_pval[idx_max];
    }
    FWER_opt = floor(idx_max+1)/n_perm;
}

void SignificantIntervalSearchWy::decrease_threshold(){
    super::decrease_threshold();
    // Recompute FWER from scratch
    int false_positives = 0; //Number of false positives (a false positive occurs if min_pval[j] <= delta)
    for(int j=0; j<n_perm; j++) false_positives += (min_pval[j]<=delta) ? 1 : 0;
    FWER = ((double)false_positives)/n_perm;
}

} /* namespace SignificantPattern */
