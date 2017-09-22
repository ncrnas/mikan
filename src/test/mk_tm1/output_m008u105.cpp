#include <iostream>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

// Overlap - 7mer-m8: 439 & 8mer: 438

class OM008U105 : public TestIOBase {
protected:
    OM008U105() {
        IFNAME1 = (char *) "mir_008.fasta";
        IFNAME2 = (char *) "utr3_105.fasta";
        O1FNAME1 = (char *) "test_output1_site_15.txt";
        O1FNAME2 = (char *) "test_output1_mrna_15.txt";
        O2FNAME1 = (char *) "test_output2_site_15.txt";
        O2FNAME2 = (char *) "test_output2_mrna_15.txt";
        OMPATH = (char *) "mk_tm1/";
    }
};

TEST_F(OM008U105, comp_site) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM008U105, comp_mrna) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}