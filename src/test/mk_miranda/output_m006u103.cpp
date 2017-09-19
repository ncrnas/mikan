#include <iostream>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class OM006U103 : public TestIOBase {
protected:
    OM006U103() {
        IFNAME1 = (char *) "mir_006.fasta";
        IFNAME2 = (char *) "utr3_103.fasta";
        O1FNAME1 = (char *) "test_output1_site_13.txt";
        O1FNAME2 = (char *) "test_output1_mrna_13.txt";
        O2FNAME1 = (char *) "test_output2_site_13.txt";
        O2FNAME2 = (char *) "test_output2_mrna_13.txt";
        OMPATH = (char *) "mk_miranda/";
    }
};

TEST_F(OM006U103, comp_site) {
    (void) mikan::MKCoreMain<mr3as::MR3Options, mr3as::MR3Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM006U103, comp_mrna) {
    (void) mikan::MKCoreMain<mr3as::MR3Options, mr3as::MR3Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}