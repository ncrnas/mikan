#include <iostream>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class SiteOut05BT2A : public TestIORH2 {
protected:
    SiteOut05BT2A() {
        IFNAME1 = (char *) "mir_004.fasta";
        IFNAME2 = (char *) "ts_05_bt_2.fasta";
        O1FNAME1 = (char *) "so05_bt_2_site_orig.txt";
        O1FNAME2 = (char *) "so05_bt_2_mrna_orig.txt";
        O2FNAME1 = (char *) "so05_bt_2_site_mk.txt";
        O2FNAME2 = (char *) "so05_bt_2_mrna_mk.txt";
        OMPATH = (char *) "mk_rh2/";
    }
};

TEST_F(SiteOut05BT2A, comp_orig_mk) {
    (void) rh2mfe::RH2CoreMain(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
    gtest_compare_two_files(o1file2, o2file2);
}

class SiteOut05BT2B : public TestIORH2 {
protected:
    SiteOut05BT2B() {
        IFNAME1 = (char *) "mir_004.fasta";
        IFNAME2 = (char *) "ts_05_bt_2.fasta";
        O1FNAME1 = (char *) "so05_bt_2_site_orig_b.txt";
        O1FNAME2 = (char *) "so05_bt_2_mrna_orig_b.txt";
        O2FNAME1 = (char *) "so05_bt_2_site_mk_b.txt";
        O2FNAME2 = (char *) "so05_bt_2_mrna_mk_b.txt";
        OMPATH = (char *) "mk_rh2/";
    }
};

TEST_F(SiteOut05BT2B, comp_orig_mk) {
    argc = 6;
    argv[5] = (char *) "--overlap=orig";

    (void) rh2mfe::RH2CoreMain(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
    gtest_compare_two_files(o1file2, o2file2);
}
}