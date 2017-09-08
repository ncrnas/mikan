#include <iostream>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class OM001U005 : public TestIOBase {
protected:
    OM001U005() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_005.fasta";
        O1FNAME1 = (char *) "test_output1_site_5.txt";
        O1FNAME2 = (char *) "test_output1_mrna_5.txt";
        O2FNAME1 = (char *) "test_output2_site_5.txt";
        O2FNAME2 = (char *) "test_output2_mrna_5.txt";
        OMPATH = (char *) "mk_tm1/";
    }
};

TEST_F(OM001U005, comp_site) {
    (void) mikan::MKCoreMain<tm1p::TM1Options, tm1p::TM1Core >(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM001U005, comp_mrna) {
    (void) mikan::MKCoreMain<tm1p::TM1Options, tm1p::TM1Core >(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}