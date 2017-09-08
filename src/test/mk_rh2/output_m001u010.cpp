#include <iostream>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

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
        OMPATH = (char *) "mk_rh2/";
    }
};

TEST_F(OM001U010, comp_site) {
    (void) mikan::MKCoreMain<rh2mfe::RH2Options, rh2mfe::RH2Core >(argc, (const char **) argv);
    gtest_compare_two_files2(o1file1, o2file1, 6, 1000, 1);
}

TEST_F(OM001U010, comp_mrna) {
    (void) mikan::MKCoreMain<rh2mfe::RH2Options, rh2mfe::RH2Core >(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}