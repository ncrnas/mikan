#include <iostream>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class SiteOut06BM1 : public TestIOBase {
protected:
    SiteOut06BM1() {
        IFNAME1 = (char *) "mir_003.fasta";
        IFNAME2 = (char *) "ts_06_bm_1.fasta";
        O1FNAME1 = (char *) "so06_bm_1_site_orig.txt";
        O1FNAME2 = (char *) "so06_bm_1_mrna_orig.txt";
        O2FNAME1 = (char *) "so06_bm_1_site_mk.txt";
        O2FNAME2 = (char *) "so06_bm_1_mrna_mk.txt";
        OMPATH = (char *) "mk_pita/";
    }
};

TEST_F(SiteOut06BM1, comp_orig_mk) {
    (void) mikan::MKCoreMain<ptddg::PITAOptions, ptddg::PITACore >(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
    gtest_compare_two_files(o1file2, o2file2);
}
}