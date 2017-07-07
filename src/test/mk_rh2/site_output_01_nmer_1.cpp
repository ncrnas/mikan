#include <iostream>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class SiteOut01Nmer1 : public TestIORH2 {
protected:
    SiteOut01Nmer1() {
        IFNAME1 = (char *) "mir_003.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";
        O1FNAME1 = (char *) "so01_nmer_1_site_orig.txt";
        O1FNAME2 = (char *) "so01_nmer_1_mrna_orig.txt";
        O2FNAME1 = (char *) "so01_nmer_1_site_mk.txt";
        O2FNAME2 = (char *) "so01_nmer_1_mrna_mk.txt";
        OMPATH = (char *) "mk_rh2/";
    }
};

TEST_F(SiteOut01Nmer1, comp_orig_mk) {
    (void) rh2mfe::RH2CoreMain(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
    gtest_compare_two_files(o1file2, o2file2);
}
}