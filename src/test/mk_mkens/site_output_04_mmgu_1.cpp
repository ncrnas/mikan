#include <iostream>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class SiteOut04MMGU1 : public TestIOBase {
protected:
    SiteOut04MMGU1() {
        IFNAME1 = (char *) "mir_003.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_1.fasta";
        O1FNAME1 = (char *) "so04_mmgu_1_site_orig.txt";
        O1FNAME2 = (char *) "so04_mmgu_1_mrna_orig.txt";
        O2FNAME1 = (char *) "so04_mmgu_1_site_mk.txt";
        O2FNAME2 = (char *) "so04_mmgu_1_mrna_mk.txt";
        OMPATH = (char *) "mk_mkens/";
    }
};

TEST_F(SiteOut04MMGU1, comp_orig_mk) {
    (void) mikan::MKCoreMain<mkens::MKEOptions, mkens::MKECore >(argc, (const char **) argv);
    gtest_compare_two_files4(o1file1, o2file1, 4, 100, 1, 5, 100, 1);
    gtest_compare_two_files2(o1file2, o2file2, 2, 100, 1);
}
}