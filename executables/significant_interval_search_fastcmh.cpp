//============================================================================
// Name        : significant_interval_search_fastcmh.cpp
// Author      : ETH
// Version     :
// Copyright   :
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include "SignificantIntervalSearchFastCmh.h"
#include <Exception.h>
#include <FeatureSet.h>
#include <Summary.h>
#include <types.h>

void usage()
{
    fprintf(stderr, "ERROR: INCORRECT SYNTAX!\n");
    fprintf(stderr, "\tUSAGE: ./program_name -eth X_file Y_file cov_file alpha L_max base_filename [out_file_pvals]\n");
    fprintf(stderr, "\tOR\n");
    fprintf(stderr,"\tUSAGE: ./program_name -plink basefilename cov_file alpha L_max base_filename [out_file_pvals]\n");
    exit(1);
}

int main(int argc, char *argv[])
{

    string xfilename, yfilename, covfilename, plinkbase, basefilename, pvalfilename;
    string filetype;
    double alpha;
    longint L_max;

    // Check input
    if(argc < 6) usage();

    int arg = 1;
    filetype = argv[arg++];
    if (filetype == "-plink")
    {
        plinkbase = argv[arg++];
    }
    else if (filetype == "-eth")
    {
        xfilename = argv[arg++];
        yfilename = argv[arg++];
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
        SignificantPattern::SignificantIntervalSearchFastCmh search;
        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_fastcmh: read\n");
        #endif
        // read covariates in one call w/ data and labels
       if (filetype == "-plink")
           search.readPlinkFilesWithCovariates(plinkbase, covfilename);
       else
           search.readETHFilesWithCovariates(xfilename, yfilename, covfilename);
        // uncomment to test separate covariates read
        // bool plinkCov = false;
        // if (filetype == "-plink") {
        //     search.readPlinkFiles(plinkbase);
        //     plinkCov = true;
        // } else {
        //     search.readETHFiles(xfilename, yfilename);
        // }
        // #ifdef DEBUG
        // fprintf(stderr,
        //         "\nsignificant_interval_search_fastcmh: read covariates '%s' in "
        //         "%s format\n",
        //         covfilename.c_str(), plinkCov ? "PLINK" : "ETH");
        // #endif
        // search.readCovariatesFile(covfilename, plinkCov);

        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_fastcmh: execute\n");
        #endif
        search.execute(alpha, L_max);
        //uncomment to test double execute
        // #ifdef DEBUG
        // fprintf(stderr, "\nsignificant_interval_search_fastcmh: execute #2\n");
        // #endif
        // search.execute(alpha, L_max); //#2

        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_fastcmh: write\n");
        #endif
        if (filetype == "-eth") {
            string data_filename = basefilename + "_data.txt";
            string label_filename = basefilename + "_label.txt";
            string cov_filename = basefilename + "_cov.txt";
            search.writeETHFilesWithCovariates(data_filename, label_filename, cov_filename);
        }
        string sig_int_filename = basefilename + "_sigints.csv";
        search.getPValsSigInts().writeToFile(sig_int_filename);
        if (pvalfilename != "") // basefilename + "_p_vals_testable_ints.txt";
            search.getPValsTestableInts().writeToFile(pvalfilename);
        string summary_filename = basefilename + "_summary.txt";
        search.getSummary().writeToFile(summary_filename);
        string timing_filename = basefilename + "_timing.txt";
        search.getProfiler().writeToFile(timing_filename);
        string filtered_filename = basefilename
                + "_sigints_filtered.corrected.csv";
        search.getFilteredIntervals().writeToFile(filtered_filename);
    }
    catch (const std::exception& e)
    {
        cout << "Caught Exception: " << e.what() << endl;
        return 1;
    }


    return 0;
}
