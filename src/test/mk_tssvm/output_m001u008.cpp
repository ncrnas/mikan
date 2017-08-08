#include <iostream>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class OM001U008 : public TestIOTSSVM {
protected:
    OM001U008() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_008.fasta";
        O1FNAME1 = (char *) "test_output1_site_8.txt";
        O1FNAME2 = (char *) "test_output1_mrna_8.txt";
        O2FNAME1 = (char *) "test_output2_site_8.txt";
        O2FNAME2 = (char *) "test_output2_mrna_8.txt";
        OMPATH = (char *) "mk_tssvm/";
    }
};

TEST_F(OM001U008, comp_site) {
    (void) tssvm::TSSVMCoreMain(argc, (const char **) argv);
    gtest_compare_two_files2(o1file1, o2file1, 5, 100, 1);
}

TEST_F(OM001U008, comp_mrna) {
    (void) tssvm::TSSVMCoreMain(argc, (const char **) argv);
    gtest_compare_two_files2(o1file2, o2file2, 2, 100, 1);
}
}