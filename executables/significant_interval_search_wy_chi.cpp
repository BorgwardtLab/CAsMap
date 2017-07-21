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

#include <SignificantIntervalSearchWyChi.h>
#include <Exception.h>
#include <FeatureSet.h>
#include <Summary.h>
#include <types.h>

void usage()
{
    fprintf(stderr, "ERROR: INCORRECT SYNTAX!\n");
    fprintf(stderr, "\tUSAGE: ./program_name [-seed n] -eth X_file Y_file alpha n_perm L_max base_filename [out_file_pvals]\n");
    fprintf(stderr, "\tOR\n");
    fprintf(stderr,"\tUSAGE: ./program_name [-seed n] -plink basefilename alpha n_perm L_max base_filename [out_file_pvals]\n");
    exit(1);

}

int main(int argc, char *argv[])
{
    string xfilename, yfilename, plinkbase, basefilename, pvalfilename;
    string filetype;
    double alpha;
    longint L_max;
    longint n_perm;

    // Check input
    if(argc < 6) usage();

    int idx = 1;
    long seed = -1;
    while (idx < argc && argv[idx][0] == '-')
    {
        if (string(argv[idx]) == "-seed")
        {
            idx++;
            seed = atol(argv[idx++]);
        }
        else if (string(argv[idx]) == "-plink" || string(argv[idx]) == "-eth")
        {
            filetype = argv[idx];
            idx++;
        }
        else usage();
    }
    if (filetype == "-eth")
    {
        xfilename = argv[idx++];
        yfilename = argv[idx++];
    }
    else if (filetype == "-plink")
        plinkbase = argv[idx++];
    else usage();
    alpha = atof(argv[idx++]);
    n_perm = atoll(argv[idx++]);
    if (n_perm < 4) {
        printf("n_perm >= 4 is required\n");
        exit(1);
    }
    L_max = atoll(argv[idx++]);
    basefilename = argv[idx++];
    if (argc > idx)
        pvalfilename = argv[idx++];

    try
    {
        SignificantPattern::SignificantIntervalSearchWyChi search;
        search.setNPerm(n_perm);
        if (seed > -1)
            search.setSeed((unsigned)seed);
        #ifdef DEBUG
        cerr << endl << "significant_itemset_search_wy_chi: read" << endl;
        #endif
        if (filetype == "-plink")
            search.readPlinkFiles(plinkbase);
        else
            search.readETHFiles(xfilename, yfilename);

        #ifdef DEBUG
        cerr << endl << "significant_itemset_search_wy_chi: execute" << endl;
        #endif
        search.execute(alpha, L_max);
        //uncomment to test:
//        cerr << "[DEBUG] RSS post-execute #1: " << search.getProfiler().getCurrMemory()/1024 << endl;
//        // * an execute error, or
//        search.setNPerm(-1);
//        try { search.execute(alpha, L_max); }
//        catch (const std::exception& e) {
//            cerr << "[DEBUG] caught Exception:\n\t" << e.what() << endl;
//        }
//        cerr << "[DEBUG] RSS post-error execute: " << search.getProfiler().getCurrMemory()/1024 << endl;
//        search.setNPerm(n_perm);
//        //* multiple executes
//        search.execute(alpha, L_max); //#2
//        cerr << "[DEBUG] RSS current post-execute #2: " << search.getProfiler().getCurrMemory()/1024 << endl;
//        search.execute(alpha, L_max); //#3
//        cerr << "[DEBUG] RSS current post-execute #3: " << search.getProfiler().getCurrMemory()/1024 << endl;

        #ifdef DEBUG
        cerr << endl << "significant_itemset_search_wy_chi: write" << endl;
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
    }
    catch (const std::exception& e)
    {
        cout << "Caught Exception: " << e.what() << endl;
        return 1;
    }


    return 0;
}
