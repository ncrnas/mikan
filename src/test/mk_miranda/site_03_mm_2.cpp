#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class Site03MM2 : public TestSiteMR3AS {
protected:
    Site03MM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_03_mm_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_miranda/";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "0";
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "0";
    }

    typedef mr3as::MR3Core<mr3as::TRNATYPE>::TIndexQGram TIdx;
    typedef mr3as::MR3Core<mr3as::TRNATYPE>::TFinder TFin;
    typedef mr3as::MR3SeedSites<mr3as::TRNATYPE> TSit;

};

TEST_F(Site03MM2, mir1_mm7) {
    read_files(false);
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(40u, sites.get_length());

    test_sites(sites, 0, "7mer_MM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_MM", 1, 24, true, -1);
    test_sites(sites, 2, "8mer_MM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MM", 3, 24, true, -1);

    test_sites(sites, 4, "MM", 4, 24, false, 0);
    test_sites(sites, 5, "MM", 5, 24, false, 0);
    test_sites(sites, 6, "7mer_MM", 6, 24, true, 0);
    test_sites(sites, 7, "7mer_MM", 7, 24, true, 0);
    test_sites(sites, 8, "8mer_MM", 8, 24, true, 0);
    test_sites(sites, 9, "8mer_MM", 9, 24, true, 0);

    test_sites(sites, 10, "MM", 10, 24, false, 0);
    test_sites(sites, 11, "MM", 11, 24, false, 0);
    test_sites(sites, 12, "7mer_MM", 12, 24, true, 1);
    test_sites(sites, 13, "7mer_MM", 13, 24, true, 1);
    test_sites(sites, 14, "8mer_MM", 14, 24, true, 1);
    test_sites(sites, 15, "8mer_MM", 15, 24, true, 1);

    test_sites(sites, 16, "MM", 16, 24, false, 0);
    test_sites(sites, 17, "MM", 17, 24, false, 0);
    test_sites(sites, 18, "7mer_MM", 18, 24, true, 2);
    test_sites(sites, 19, "7mer_MM", 19, 24, true, 2);
    test_sites(sites, 20, "8mer_MM", 20, 24, true, 2);
    test_sites(sites, 21, "8mer_MM", 21, 24, true, 2);

    test_sites(sites, 22, "MM", 22, 24, false, 0);
    test_sites(sites, 23, "MM", 23, 24, false, 0);
    test_sites(sites, 24, "7mer_MM", 24, 24, true, 3);
    test_sites(sites, 25, "7mer_MM", 25, 24, true, 3);
    test_sites(sites, 26, "8mer_MM", 26, 24, true, 3);
    test_sites(sites, 27, "8mer_MM", 27, 24, true, 3);

    test_sites(sites, 28, "MM", 28, 24, false, 0);
    test_sites(sites, 29, "MM", 29, 24, false, 0);
    test_sites(sites, 30, "7mer_MM", 30, 24, true, 4);
    test_sites(sites, 31, "7mer_MM", 31, 24, true, 4);
    test_sites(sites, 32, "8mer_MM", 32, 24, true, 4);
    test_sites(sites, 33, "8mer_MM", 33, 24, true, 4);

    test_sites(sites, 34, "MM", 34, 24, false, 0);
    test_sites(sites, 35, "MM", 35, 24, false, 0);
    test_sites(sites, 36, "7mer_MM", 36, 24, true, 5);
    test_sites(sites, 37, "7mer_MM", 37, 24, true, 5);
    test_sites(sites, 38, "8mer_MM", 38, 24, true, 5);
    test_sites(sites, 39, "8mer_MM", 39, 24, true, 5);
}

TEST_F(Site03MM2, mir1_mm8) {
    read_files(false);
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[4] = "0:1";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(40u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "8mer_MM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MM", 3, 24, true, -1);

    test_sites(sites, 4, "MM", 4, 24, false, 0);
    test_sites(sites, 5, "MM", 5, 24, false, 0);
    test_sites(sites, 6, "MM", 6, 24, false, 0);
    test_sites(sites, 7, "MM", 7, 24, false, 0);
    test_sites(sites, 8, "8mer_MM", 8, 24, true, 0);
    test_sites(sites, 9, "8mer_MM", 9, 24, true, 0);

    test_sites(sites, 10, "MM", 10, 24, false, 0);
    test_sites(sites, 11, "MM", 11, 24, false, 0);
    test_sites(sites, 12, "MM", 12, 24, false, 0);
    test_sites(sites, 13, "MM", 13, 24, false, 0);
    test_sites(sites, 14, "8mer_MM", 14, 24, true, 1);
    test_sites(sites, 15, "8mer_MM", 15, 24, true, 1);

    test_sites(sites, 16, "MM", 16, 24, false, 0);
    test_sites(sites, 17, "MM", 17, 24, false, 0);
    test_sites(sites, 18, "MM", 18, 24, false, 0);
    test_sites(sites, 19, "MM", 19, 24, false, 0);
    test_sites(sites, 20, "8mer_MM", 20, 24, true, 2);
    test_sites(sites, 21, "8mer_MM", 21, 24, true, 2);

    test_sites(sites, 22, "MM", 22, 24, false, 0);
    test_sites(sites, 23, "MM", 23, 24, false, 0);
    test_sites(sites, 24, "MM", 24, 24, false, 0);
    test_sites(sites, 25, "MM", 25, 24, false, 0);
    test_sites(sites, 26, "8mer_MM", 26, 24, true, 3);
    test_sites(sites, 27, "8mer_MM", 27, 24, true, 3);

    test_sites(sites, 28, "MM", 28, 24, false, 0);
    test_sites(sites, 29, "MM", 29, 24, false, 0);
    test_sites(sites, 30, "MM", 30, 24, false, 0);
    test_sites(sites, 31, "MM", 31, 24, false, 0);
    test_sites(sites, 32, "8mer_MM", 32, 24, true, 4);
    test_sites(sites, 33, "8mer_MM", 33, 24, true, 4);

    test_sites(sites, 34, "MM", 34, 24, false, 0);
    test_sites(sites, 35, "MM", 35, 24, false, 0);
    test_sites(sites, 36, "MM", 36, 24, false, 0);
    test_sites(sites, 37, "MM", 37, 24, false, 0);
    test_sites(sites, 38, "8mer_MM", 38, 24, true, 5);
    test_sites(sites, 39, "8mer_MM", 39, 24, true, 5);
}

TEST_F(Site03MM2, mir1_def) {
    read_files(false);
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "+";
    mSeedDef[4] = "1:1";
    mSeedDef[5] = "1";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(106u, sites.get_length());

    test_sites(sites, 0, "7mer_MM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_MM", 1, 24, true, -1);
    test_sites(sites, 2, "8mer_MM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MM", 3, 24, true, -1);

    test_sites(sites, 4, "MM", 4, 24, false, 0);
    test_sites(sites, 5, "MM", 5, 24, false, 0);
    test_sites(sites, 6, "7mer_MM", 6, 24, true, 0);
    test_sites(sites, 7, "7mer_MM", 7, 24, true, 0);
    test_sites(sites, 8, "8mer_MM", 8, 24, true, 0);
    test_sites(sites, 9, "8mer_MM", 9, 24, true, 0);

    test_sites(sites, 10, "MM", 10, 24, false, 0);
    test_sites(sites, 11, "MM", 11, 24, false, 0);
    test_sites(sites, 12, "7mer_MM", 12, 24, true, 1);
    test_sites(sites, 13, "7mer_MM", 13, 24, true, 1);
    test_sites(sites, 14, "8mer_MM", 14, 24, true, 1);
    test_sites(sites, 15, "8mer_MM", 15, 24, true, 1);

    test_sites(sites, 16, "MM", 16, 24, false, 0);
    test_sites(sites, 17, "MM", 17, 24, false, 0);
    test_sites(sites, 18, "7mer_MM", 18, 24, true, 2);
    test_sites(sites, 19, "7mer_MM", 19, 24, true, 2);
    test_sites(sites, 20, "8mer_MM", 20, 24, true, 2);
    test_sites(sites, 21, "8mer_MM", 21, 24, true, 2);

    test_sites(sites, 22, "MM", 22, 24, false, 0);
    test_sites(sites, 23, "MM", 23, 24, false, 0);
    test_sites(sites, 24, "7mer_MM", 24, 24, true, 3);
    test_sites(sites, 25, "7mer_MM", 25, 24, true, 3);
    test_sites(sites, 26, "8mer_MM", 26, 24, true, 3);
    test_sites(sites, 27, "8mer_MM", 27, 24, true, 3);

    test_sites(sites, 28, "MM", 28, 24, false, 0);
    test_sites(sites, 29, "MM", 29, 24, false, 0);
    test_sites(sites, 30, "7mer_MM", 30, 24, true, 4);
    test_sites(sites, 31, "7mer_MM", 31, 24, true, 4);
    test_sites(sites, 32, "8mer_MM", 32, 24, true, 4);
    test_sites(sites, 33, "8mer_MM", 33, 24, true, 4);

    test_sites(sites, 34, "MM", 34, 24, false, 0);
    test_sites(sites, 35, "MM", 35, 24, false, 0);
    test_sites(sites, 36, "7mer_MM", 36, 24, true, 5);
    test_sites(sites, 37, "7mer_MM", 37, 24, true, 5);
    test_sites(sites, 38, "8mer_MM", 38, 24, true, 5);
    test_sites(sites, 39, "8mer_MM", 39, 24, true, 5);

    test_sites(sites, 40, "7mer_MMGU", 4, 23, true, 2);
    test_sites(sites, 41, "7mer_MMGU", 5, 23, true, 2);

    test_sites(sites, 42, "7mer_MMGU", 0, 12, true, 3);
    test_sites(sites, 43, "7mer_MMGU", 1, 12, true, 3);
    test_sites(sites, 44, "7mer_MMGU", 2, 12, true, 3);
    test_sites(sites, 45, "7mer_MMGU", 3, 12, true, 3);

    test_sites(sites, 46, "7mer_MMGU", 4, 12, true, 3);
    test_sites(sites, 47, "7mer_MMGU", 5, 12, true, 3);
    test_sites(sites, 48, "7mer_MMGU", 6, 12, true, 3);
    test_sites(sites, 49, "7mer_MMGU", 7, 12, true, 3);

    test_sites(sites, 50, "7mer_MMGU", 8, 12, true, 3);
    test_sites(sites, 51, "7mer_MMGU", 9, 12, true, 3);
    test_sites(sites, 52, "7mer_MMGU", 10, 12, true, 3);
    test_sites(sites, 53, "7mer_MMGU", 11, 12, true, 3);
    test_sites(sites, 54, "7mer_MMGU", 12, 12, true, 3);

    test_sites(sites, 55, "7mer_MMGU", 13, 12, true, 3);
    test_sites(sites, 56, "7mer_MMGU", 14, 12, true, 3);
    test_sites(sites, 57, "7mer_MMGU", 15, 12, true, 3);
    test_sites(sites, 58, "7mer_MMGU", 16, 12, true, 3);

    test_sites(sites, 59, "7mer_MMGU", 17, 12, true, 3);
    test_sites(sites, 60, "7mer_MMGU", 18, 12, true, 3);
    test_sites(sites, 61, "7mer_MMGU", 19, 12, true, 3);
    test_sites(sites, 62, "7mer_MMGU", 20, 12, true, 3);
    test_sites(sites, 63, "7mer_MMGU", 21, 12, true, 3);

    test_sites(sites, 64, "7mer_MMGU", 22, 12, true, 3);
    test_sites(sites, 65, "7mer_MMGU", 23, 12, true, 3);
    test_sites(sites, 66, "7mer_MMGU", 24, 12, true, 3);
    test_sites(sites, 67, "7mer_MMGU", 25, 12, true, 3);

    test_sites(sites, 68, "7mer_MMGU", 26, 12, true, 3);
    test_sites(sites, 69, "7mer_MMGU", 27, 12, true, 3);
    test_sites(sites, 70, "7mer_MMGU", 28, 12, true, 3);
    test_sites(sites, 71, "7mer_MMGU", 29, 12, true, 3);
    test_sites(sites, 72, "7mer_MMGU", 30, 12, true, 3);

    test_sites(sites, 73, "7mer_MMGU", 31, 12, true, 3);
    test_sites(sites, 74, "7mer_MMGU", 32, 12, true, 3);
    test_sites(sites, 75, "7mer_MMGU", 33, 12, true, 3);
    test_sites(sites, 76, "7mer_MMGU", 34, 12, true, 3);

    test_sites(sites, 77, "7mer_MMGU", 35, 12, true, 3);
    test_sites(sites, 78, "7mer_MMGU", 36, 12, true, 3);
    test_sites(sites, 79, "7mer_MMGU", 37, 12, true, 3);
    test_sites(sites, 80, "7mer_MMGU", 38, 12, true, 3);
    test_sites(sites, 81, "7mer_MMGU", 39, 12, true, 3);

    test_sites(sites, 82, "BT", 36, 23, false, 0);
    test_sites(sites, 83, "BT", 37, 23, false, 0);
    test_sites(sites, 84, "BT", 38, 23, false, 0);
    test_sites(sites, 85, "BT", 39, 23, false, 0);
    test_sites(sites, 86, "7mer_BT", 0, 23, true, 0);
    test_sites(sites, 87, "7mer_BT", 1, 23, true, 0);
    test_sites(sites, 88, "BT", 2, 23, false, 0);
    test_sites(sites, 89, "BT", 3, 23, false, 0);
    test_sites(sites, 90, "BT", 34, 23, false, 0);
    test_sites(sites, 91, "BT", 35, 23, false, 0);
    test_sites(sites, 92, "7mer_BT", 4, 23, true, 1);
    test_sites(sites, 93, "7mer_BT", 5, 23, true, 1);

    test_sites(sites, 94, "BT", 28, 24, false, 0);
    test_sites(sites, 95, "BT", 29, 24, false, 0);
    test_sites(sites, 96, "7mer_BT", 30, 24, true, 4);
    test_sites(sites, 97, "BT", 31, 24, false, 0);
    test_sites(sites, 98, "8mer_BT", 32, 24, true, 4);
    test_sites(sites, 99, "BT", 33, 24, false, 0);
    test_sites(sites, 100, "BT", 34, 24, false, 0);
    test_sites(sites, 101, "BT", 35, 24, false, 0);
    test_sites(sites, 102, "7mer_BT", 36, 24, true, 5);
    test_sites(sites, 103, "BT", 37, 24, false, 0);
    test_sites(sites, 104, "8mer_BT", 38, 24, true, 5);
    test_sites(sites, 105, "BT", 39, 24, false, 0);
}
}
