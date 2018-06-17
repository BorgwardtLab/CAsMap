//============================================================================
// Name        : significant_itemset_search_facs.cpp
// Author      : ETH
// Version     :
// Copyright   :
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include "SignificantItemsetSearchFacs.h"
#include <Exception.h>
#include <FeatureSet.h>
#include <Summary.h>
#include <types.h>

void usage()
{
    fprintf(stderr, "ERROR: INCORRECT SYNTAX!\n");
    fprintf(stderr, "\nTO USE DOMINANT ENCODING\n");
    fprintf(stderr, "\tUSAGE: ./program_name -eth X_file Y_file cov_file alpha L_max base_filename [out_file_pvals]\n");
    fprintf(stderr, "\tOR\n");
    fprintf(stderr,"\tUSAGE: ./program_name -plink basefilename cov_file alpha L_max base_filename [out_file_pvals]\n");
    fprintf(stderr, "\nTO USE RECESSIVE ENCODING\n");
    fprintf(stderr, "\tUSAGE: ./program_name -eth_recessive cov_file X_file Y_file alpha L_max base_filename [out_file_pvals]\n");
    fprintf(stderr, "\tOR\n");
    fprintf(stderr,"\tUSAGE: ./program_name -plink_recessive cov_file basefilename alpha L_max base_filename [out_file_pvals]\n");
    exit(1);
}

int main(int argc, char *argv[])
{

    string xfilename, yfilename, covfilename, plinkbase, basefilename, pvalfilename;
    string filetype;
    string encoding;
    double alpha;
    longint L_max;

    // Check input
    if(argc < 6) usage();

    int arg = 1;
    filetype = argv[arg++];
    if (filetype == "-plink" or filetype == "-plink_dominant")
    {
        plinkbase = argv[arg++];
        encoding = "dominant";
    }
    else if (filetype == "-plink_recesssive"){
        plinkbase = argv[arg++];
        encoding = "recessive";
    }
    else if (filetype == "-eth" or filetype == "-eth_dominant")
    {
        xfilename = argv[arg++];
        yfilename = argv[arg++];
        encoding = "dominant";
    }
    else if (filetype == "-eth_recessive"){
        xfilename = argv[arg++];
        yfilename = argv[arg++];
        encoding = "recessive";
    }
    else usage();
    covfilename = argv[arg++];
    alpha = atof(argv[arg++]);
    L_max = atoll(argv[arg++]);
    basefilename = argv[arg++];

    if (arg < argc)
        pvalfilename = argv[arg++];

    try
    {
        SignificantPattern::SignificantItemsetSearchFacs search;
        //uncomment to test covariates read on one call
        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_itemset_search_facs: read\n");
        #endif
//        if (filetype == "-plink")
//            search.readPlinkFilesWithCovariates(plinkbase, covfilename);
//        else
//            search.readETHFilesWithCovariates(xfilename, yfilename, covfilename);
        //comment out covariates read to test a default, single covariate
        if (filetype == "-plink")
            search.readPlinkFiles(plinkbase, encoding);
        else
            search.readETHFiles(xfilename, yfilename, encoding);
        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_itemset_search_facs: read covariates '%s'\n", covfilename.c_str());
        #endif
        search.readCovariatesFile(covfilename);

        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_itemset_search_facs: execute\n");
        #endif
        search.execute(alpha, L_max);
        //uncomment to test double execute
        // #ifdef DEBUG
        // fprintf(stderr, "significant_itemset_search_facs: execute #2\n");
        // #endif
        // search.execute(alpha, L_max); //#2

        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_itemset_search_facs: write\n");
        #endif
        if (filetype == "-eth") {
            string data_filename = basefilename + "_data.txt";
            string label_filename = basefilename + "_label.txt";
            string cov_filename = basefilename + "_cov.txt";
            search.writeETHFilesWithCovariates(data_filename, label_filename, cov_filename);
        }
        string sig_int_filename = basefilename + "_sigisets.csv";
        search.getPValsSigIsets().writeToFile(sig_int_filename);
        if (pvalfilename != "") // basefilename + "_p_vals_testable_isets.txt";
            search.getPValsTestableIsets().writeToFile(pvalfilename);
        string summary_filename = basefilename + "_summary.txt";
        search.getSummary().writeToFile(summary_filename);
        string timing_filename = basefilename + "_timing.txt";
        search.getProfiler().writeToFile(timing_filename);
    }
    catch (const std::exception& e)
    {
        cout << "Caught Exception: " << e.what() << endl;
        return 1;
    }


    return 0;
}
