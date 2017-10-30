#include <iostream>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

// RNA5 with Ns

class OM005U102 : public TestIOBase {
protected:
    OM005U102() {
        IFNAME1 = (char *) "mir_005.fasta";
        IFNAME2 = (char *) "utr3_102.fasta";
        O1FNAME1 = (char *) "test_output1_site_12.txt";
        O1FNAME2 = (char *) "test_output1_mrna_12.txt";
        O2FNAME1 = (char *) "test_output2_site_12.txt";
        O2FNAME2 = (char *) "test_output2_mrna_12.txt";
        OMPATH = (char *) "mk_pita/";
    }
};

TEST_F(OM005U102, comp_site) {
    (void) mikan::MKMain<ptddg::PITAOptions, ptddg::PITACore>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM005U102, comp_mrna) {
    (void) mikan::MKMain<ptddg::PITAOptions, ptddg::PITACore>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}