#include <iostream>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

    class SiteOut05BT1 : public TestIOTM1
    {
    protected:
        SiteOut05BT1() {
            IFNAME1 = (char *)"mir_003.fasta";
            IFNAME2 = (char *)"ts_05_bt_1.fasta";
            O1FNAME1 = (char *)"so05_bt_1_site_orig.txt";
            O1FNAME2 = (char *)"so05_bt_1_mrna_orig.txt";
            O2FNAME1 = (char *)"so05_bt_1_site_mk.txt";
            O2FNAME2 = (char *)"so05_bt_1_mrna_mk.txt";
            OMPATH = (char *)"mk_tm1/";
        }
    };

    TEST_F(SiteOut05BT1, comp_orig_mk) {
        (void)tm1p::TM1CoreMain(argc, (const char **)argv);
        gtest_compare_two_files(o1file1, o2file1);
        // TODO Need to add sort, index, and round functions
//        gtest_compare_two_files(o1file2, o2file2);
    }
}