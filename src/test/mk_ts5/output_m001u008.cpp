#include <iostream>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class OM001U008 : public TestIOBase {
protected:
    OM001U008() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_008.fasta";
        O1FNAME1 = (char *) "test_output1_site_8.txt";
        O1FNAME2 = (char *) "test_output1_mrna_8.txt";
        O2FNAME1 = (char *) "test_output2_site_8.txt";
        O2FNAME2 = (char *) "test_output2_mrna_8.txt";
        OMPATH = (char *) "mk_ts5/";
    }
};

TEST_F(OM001U008, comp_site) {
    (void) mikan::MKCoreMain<ts5cs::TS5Options, ts5cs::TS5Core >(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM001U008, comp_mrna) {
    (void) mikan::MKCoreMain<ts5cs::TS5Options, ts5cs::TS5Core >(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}