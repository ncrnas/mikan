#include <iostream>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

// NA feature values and score:-10000

class OM010U109 : public TestIOBase {
protected:
    OM010U109() {
        IFNAME1 = (char *) "mir_010.fasta";
        IFNAME2 = (char *) "utr3_109.fasta";
        O1FNAME1 = (char *) "test_output1_site_19.txt";
        O1FNAME2 = (char *) "test_output1_mrna_19.txt";
        O2FNAME1 = (char *) "test_output2_site_19.txt";
        O2FNAME2 = (char *) "test_output2_mrna_19.txt";
        OMPATH = (char *) "mk_tm1/";
    }
};

TEST_F(OM010U109, comp_site) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM010U109, comp_mrna) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}