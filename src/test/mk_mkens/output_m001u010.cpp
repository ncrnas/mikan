#include <iostream>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class OM001U010 : public TestIOBase {
protected:
    OM001U010() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_010.fasta";
        O1FNAME1 = (char *) "test_output1_site_10.txt";
        O1FNAME2 = (char *) "test_output1_mrna_10.txt";
        O2FNAME1 = (char *) "test_output2_site_10.txt";
        O2FNAME2 = (char *) "test_output2_mrna_10.txt";
        OMPATH = (char *) "mk_mkens/";
    }
};

TEST_F(OM001U010, comp_site) {
    (void) mikan::MKCoreMain<mkens::MKEOptions, mkens::MKECore >(argc, (const char **) argv);
    gtest_compare_two_files4(o1file1, o2file1, 4, 100, 1, 5, 100, 1);
}

TEST_F(OM001U010, comp_mrna) {
    (void) mikan::MKCoreMain<mkens::MKEOptions, mkens::MKECore >(argc, (const char **) argv);
    gtest_compare_two_files2(o1file2, o2file2, 2, 100, 1);
}
}