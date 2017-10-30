#include <iostream>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class OM001U007 : public TestIOBase {
protected:
    OM001U007() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_007.fasta";
        O1FNAME1 = (char *) "test_output1_site_7.txt";
        O1FNAME2 = (char *) "test_output1_mrna_7.txt";
        O2FNAME1 = (char *) "test_output2_site_7.txt";
        O2FNAME2 = (char *) "test_output2_mrna_7.txt";
        OMPATH = (char *) "mk_miranda/";
    };
};

TEST_F(OM001U007, comp_site) {
    (void) mikan::MKMain<mr3as::MR3Options, mr3as::MR3Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM001U007, comp_mrna) {
    (void) mikan::MKMain<mr3as::MR3Options, mr3as::MR3Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}