#include <iostream>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

    class OM001U007 : public TestIOTSSVM
    {
    protected:
        OM001U007() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_007.fasta";
            O1FNAME1 = (char *)"test_output1_site_7.txt";
            O1FNAME2 = (char *)"test_output1_mrna_7.txt";
            O2FNAME1 = (char *)"test_output2_site_7.txt";
            O2FNAME2 = (char *)"test_output2_mrna_7.txt";
            OMPATH = (char *)"mk_tssvm/";
        };
    };

    TEST_F(OM001U007, comp_site) {
        (void)tssvm::TSSVMCoreMain(argc, (const char **)argv);
        gtest_compare_two_files(o1file1, o2file1);
    }

    TEST_F(OM001U007, comp_mrna) {
        (void)tssvm::TSSVMCoreMain(argc, (const char **)argv);
        gtest_compare_two_files(o1file2, o2file2);
    }
}