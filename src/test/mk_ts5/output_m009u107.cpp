#include <iostream>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

// Overlap - 8mer: 1551 & 7mer-A1: 551

class OM009U107 : public TestIOBase {
protected:
    OM009U107() {
        IFNAME1 = (char *) "mir_009.fasta";
        IFNAME2 = (char *) "utr3_107.fasta";
        O1FNAME1 = (char *) "test_output1_site_17.txt";
        O1FNAME2 = (char *) "test_output1_mrna_17.txt";
        O2FNAME1 = (char *) "test_output2_site_17.txt";
        O2FNAME2 = (char *) "test_output2_mrna_17.txt";
        OMPATH = (char *) "mk_ts5/";
    }
};

TEST_F(OM009U107, comp_site) {
    (void) mikan::MKMain<ts5cs::TS5Options, ts5cs::TS5Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM009U107, comp_mrna) {
    (void) mikan::MKMain<ts5cs::TS5Options, ts5cs::TS5Core>(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}