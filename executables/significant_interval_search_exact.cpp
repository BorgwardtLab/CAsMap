//============================================================================
// Name        : significant_interval_search_exact.cpp
// Author      : ETH
// Version     :
// Copyright   :
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include <SignificantIntervalSearchExact.h>
#include <Exception.h>
#include <FeatureSet.h>
#include <Summary.h>
#include <types.h>

void usage()
{
    fprintf(stderr, "ERROR: INCORRECT SYNTAX!\n");
    fprintf(stderr, "\tUSAGE: ./program_name -eth X_file Y_file alpha L_max base_filename [out_file_pvals]\n");
    fprintf(stderr, "\tOR\n");
    fprintf(stderr,"\tUSAGE: ./program_name -plink basefilename alpha L_max base_filename [out_file_pvals]\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    string xfilename, yfilename, plinkbase, basefilename, pvalfilename;
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
    alpha = atof(argv[arg++]);
    L_max = atoll(argv[arg++]);
    basefilename = argv[arg++];

    if (arg < argc)
        pvalfilename = argv[arg++];

    try
    {
        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_exact: init\n");
        #endif
        SignificantPattern::SignificantIntervalSearchExact search;
        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_exact: read\n");
        #endif
        if (filetype == "-plink")
            search.readPlinkFiles(plinkbase);
        else
            search.readETHFiles(xfilename, yfilename);

        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_exact: execute\n");
        #endif
        search.execute(alpha, L_max);
        //uncomment to test double execute
        // #ifdef DEBUG
        // fprintf(stderr, "\nsignificant_interval_search_exact: execute #2\n");
        // #endif
        // search.execute(alpha, L_max); //#2

        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_exact: write\n");
        #endif
        if (filetype == "-eth") {
            string data_filename = basefilename + "_data.txt";
            string label_filename = basefilename + "_label.txt";
            search.writeETHFiles(data_filename, label_filename);
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
        #ifdef DEBUG
        fprintf(stderr, "\nsignificant_interval_search_exact: clean\n");
        #endif
    }
    catch (const std::exception& e)
    {
        cout << "Caught Exception: " << e.what() << endl;
        return 1;
    }


    return 0;
}
