#include <iostream>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

// Incorrect feature values

class OM003U106 : public TestIOBase {
protected:
    OM003U106() {
        IFNAME1 = (char *) "mir_003.fasta";
        IFNAME2 = (char *) "utr3_106.fasta";
        O1FNAME1 = (char *) "test_output1_site_16.txt";
        O1FNAME2 = (char *) "test_output1_mrna_16.txt";
        O2FNAME1 = (char *) "test_output2_site_16.txt";
        O2FNAME2 = (char *) "test_output2_mrna_16.txt";
        OMPATH = (char *) "mk_tm1/";
    }
};

TEST_F(OM003U106, comp_site) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM003U106, comp_mrna) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}