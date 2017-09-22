#include <iostream>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

// Incorrect upstream substring

class OM004U108 : public TestIOBase {
protected:
    OM004U108() {
        IFNAME1 = (char *) "mir_004.fasta";
        IFNAME2 = (char *) "utr3_108.fasta";
        O1FNAME1 = (char *) "test_output1_site_18.txt";
        O1FNAME2 = (char *) "test_output1_mrna_18.txt";
        O2FNAME1 = (char *) "test_output2_site_18.txt";
        O2FNAME2 = (char *) "test_output2_mrna_18.txt";
        OMPATH = (char *) "mk_ts5/";
    }
};

TEST_F(OM004U108, comp_site) {
    (void) mikan::MKCoreMain<ts5cs::TS5Options, ts5cs::TS5Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM004U108, comp_mrna) {
    (void) mikan::MKCoreMain<ts5cs::TS5Options, ts5cs::TS5Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}