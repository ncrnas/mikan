#include <iostream>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

    class OM001U009 : public TestIOTSSVM
    {
    protected:
        OM001U009() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_009.fasta";
            O1FNAME1 = (char *)"test_output1_site_9.txt";
            O1FNAME2 = (char *)"test_output1_mrna_9.txt";
            O2FNAME1 = (char *)"test_output2_site_9.txt";
            O2FNAME2 = (char *)"test_output2_mrna_9.txt";
            OMPATH = (char *)"mk_tssvm/";
        }
    };

    TEST_F(OM001U009, comp_site) {
        (void)tssvm::TSSVMCoreMain(argc, (const char **)argv);
        gtest_compare_two_files(o1file1, o2file1);
    }

    TEST_F(OM001U009, comp_mrna) {
        (void)tssvm::TSSVMCoreMain(argc, (const char **)argv);
        gtest_compare_two_files(o1file2, o2file2);
    }
}