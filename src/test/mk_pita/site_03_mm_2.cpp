#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site03MM2 : public TestSitePITA {
protected:
    Site03MM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_03_mm_2.fasta";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "0";
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "0";
    }

};

TEST_F(Site03MM2, mir1_mm7) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(28u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "8mer_MM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MM", 3, 24, true, -1);

//    test_sites(sites, 4, "MM", 4, 24, false, 0);
//    test_sites(sites, 5, "MM", 5, 24, false, 0);
    test_sites(sites, 4, "7mer_MM", 6, 24, true, 0);
    test_sites(sites, 5, "7mer_MM", 7, 24, true, 0);
    test_sites(sites, 6, "8mer_MM", 8, 24, true, 0);
    test_sites(sites, 7, "8mer_MM", 9, 24, true, 0);

//    test_sites(sites, 10, "MM", 10, 24, false, 0);
//    test_sites(sites, 11, "MM", 11, 24, false, 0);
    test_sites(sites, 8, "7mer_MM", 12, 24, true, 1);
    test_sites(sites, 9, "7mer_MM", 13, 24, true, 1);
    test_sites(sites, 10, "8mer_MM", 14, 24, true, 1);
    test_sites(sites, 11, "8mer_MM", 15, 24, true, 1);

//    test_sites(sites, 16, "MM", 16, 24, false, 0);
//    test_sites(sites, 17, "MM", 17, 24, false, 0);
    test_sites(sites, 12, "7mer_MM", 18, 24, true, 2);
    test_sites(sites, 13, "7mer_MM", 19, 24, true, 2);
    test_sites(sites, 14, "8mer_MM", 20, 24, true, 2);
    test_sites(sites, 15, "8mer_MM", 21, 24, true, 2);

//    test_sites(sites, 22, "MM", 22, 24, false, 0);
//    test_sites(sites, 23, "MM", 23, 24, false, 0);
    test_sites(sites, 16, "7mer_MM", 24, 24, true, 3);
    test_sites(sites, 17, "7mer_MM", 25, 24, true, 3);
    test_sites(sites, 18, "8mer_MM", 26, 24, true, 3);
    test_sites(sites, 19, "8mer_MM", 27, 24, true, 3);

//    test_sites(sites, 28, "MM", 28, 24, false, 0);
//    test_sites(sites, 29, "MM", 29, 24, false, 0);
    test_sites(sites, 20, "7mer_MM", 30, 24, true, 4);
    test_sites(sites, 21, "7mer_MM", 31, 24, true, 4);
    test_sites(sites, 22, "8mer_MM", 32, 24, true, 4);
    test_sites(sites, 23, "8mer_MM", 33, 24, true, 4);

//    test_sites(sites, 34, "MM", 34, 24, false, 0);
//    test_sites(sites, 35, "MM", 35, 24, false, 0);
    test_sites(sites, 24, "7mer_MM", 36, 24, true, 5);
    test_sites(sites, 25, "7mer_MM", 37, 24, true, 5);
    test_sites(sites, 26, "8mer_MM", 38, 24, true, 5);
    test_sites(sites, 27, "8mer_MM", 39, 24, true, 5);
}

TEST_F(Site03MM2, mir1_mm8) {
    mSeedDef[4] = "0:1";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(16u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "8mer_MM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MM", 3, 24, true, -1);

//    test_sites(sites, 4, "MM", 4, 24, false, 0);
//    test_sites(sites, 5, "MM", 5, 24, false, 0);
//    test_sites(sites, 6, "MM", 6, 24, false, 0);
//    test_sites(sites, 7, "MM", 7, 24, false, 0);
    test_sites(sites, 4, "8mer_MM", 8, 24, true, 0);
    test_sites(sites, 5, "8mer_MM", 9, 24, true, 0);

//    test_sites(sites, 10, "MM", 10, 24, false, 0);
//    test_sites(sites, 11, "MM", 11, 24, false, 0);
//    test_sites(sites, 12, "MM", 12, 24, false, 0);
//    test_sites(sites, 13, "MM", 13, 24, false, 0);
    test_sites(sites, 6, "8mer_MM", 14, 24, true, 1);
    test_sites(sites, 7, "8mer_MM", 15, 24, true, 1);

//    test_sites(sites, 16, "MM", 16, 24, false, 0);
//    test_sites(sites, 17, "MM", 17, 24, false, 0);
//    test_sites(sites, 18, "MM", 18, 24, false, 0);
//    test_sites(sites, 19, "MM", 19, 24, false, 0);
    test_sites(sites, 8, "8mer_MM", 20, 24, true, 2);
    test_sites(sites, 9, "8mer_MM", 21, 24, true, 2);

//    test_sites(sites, 22, "MM", 22, 24, false, 0);
//    test_sites(sites, 23, "MM", 23, 24, false, 0);
//    test_sites(sites, 24, "MM", 24, 24, false, 0);
//    test_sites(sites, 25, "MM", 25, 24, false, 0);
    test_sites(sites, 10, "8mer_MM", 26, 24, true, 3);
    test_sites(sites, 11, "8mer_MM", 27, 24, true, 3);

//    test_sites(sites, 28, "MM", 28, 24, false, 0);
//    test_sites(sites, 29, "MM", 29, 24, false, 0);
//    test_sites(sites, 30, "MM", 30, 24, false, 0);
//    test_sites(sites, 31, "MM", 31, 24, false, 0);
    test_sites(sites, 12, "8mer_MM", 32, 24, true, 4);
    test_sites(sites, 13, "8mer_MM", 33, 24, true, 4);

//    test_sites(sites, 34, "MM", 34, 24, false, 0);
//    test_sites(sites, 35, "MM", 35, 24, false, 0);
//    test_sites(sites, 36, "MM", 36, 24, false, 0);
//    test_sites(sites, 37, "MM", 37, 24, false, 0);
    test_sites(sites, 14, "8mer_MM", 38, 24, true, 5);
    test_sites(sites, 15, "8mer_MM", 39, 24, true, 5);
}

TEST_F(Site03MM2, mir1_def) {
    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(16u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "8mer_MM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MM", 3, 24, true, -1);

//    test_sites(sites, 4, "MM", 4, 24, false, 0);
//    test_sites(sites, 5, "MM", 5, 24, false, 0);
//    test_sites(sites, 6, "MM", 6, 24, false, 0);
//    test_sites(sites, 7, "MM", 7, 24, false, 0);
    test_sites(sites, 4, "8mer_MM", 8, 24, true, 0);
    test_sites(sites, 5, "8mer_MM", 9, 24, true, 0);

//    test_sites(sites, 10, "MM", 10, 24, false, 0);
//    test_sites(sites, 11, "MM", 11, 24, false, 0);
//    test_sites(sites, 12, "MM", 12, 24, false, 0);
//    test_sites(sites, 13, "MM", 13, 24, false, 0);
    test_sites(sites, 6, "8mer_MM", 14, 24, true, 1);
    test_sites(sites, 7, "8mer_MM", 15, 24, true, 1);

//    test_sites(sites, 16, "MM", 16, 24, false, 0);
//    test_sites(sites, 17, "MM", 17, 24, false, 0);
//    test_sites(sites, 18, "MM", 18, 24, false, 0);
//    test_sites(sites, 19, "MM", 19, 24, false, 0);
    test_sites(sites, 8, "8mer_MM", 20, 24, true, 2);
    test_sites(sites, 9, "8mer_MM", 21, 24, true, 2);

//    test_sites(sites, 22, "MM", 22, 24, false, 0);
//    test_sites(sites, 23, "MM", 23, 24, false, 0);
//    test_sites(sites, 24, "MM", 24, 24, false, 0);
//    test_sites(sites, 25, "MM", 25, 24, false, 0);
    test_sites(sites, 10, "8mer_MM", 26, 24, true, 3);
    test_sites(sites, 11, "8mer_MM", 27, 24, true, 3);

//    test_sites(sites, 28, "MM", 28, 24, false, 0);
//    test_sites(sites, 29, "MM", 29, 24, false, 0);
//    test_sites(sites, 30, "MM", 30, 24, false, 0);
//    test_sites(sites, 31, "MM", 31, 24, false, 0);
    test_sites(sites, 12, "8mer_MM", 32, 24, true, 4);
    test_sites(sites, 13, "8mer_MM", 33, 24, true, 4);

//    test_sites(sites, 34, "MM", 34, 24, false, 0);
//    test_sites(sites, 35, "MM", 35, 24, false, 0);
//    test_sites(sites, 36, "MM", 36, 24, false, 0);
//    test_sites(sites, 37, "MM", 37, 24, false, 0);
    test_sites(sites, 14, "8mer_MM", 38, 24, true, 5);
    test_sites(sites, 15, "8mer_MM", 39, 24, true, 5);

//    test_sites(sites, 40, "MMGU", 4, 23, false, 0);
//    test_sites(sites, 41, "MMGU", 5, 23, false, 0);

//    test_sites(sites, 42, "", 0, 12, false, 0);
//    test_sites(sites, 43, "", 1, 12, false, 0);
//    test_sites(sites, 44, "", 2, 12, false, 0);
//    test_sites(sites, 45, "", 3, 12, false, 0);
//    test_sites(sites, 46, "", 4, 12, false, 0);
//    test_sites(sites, 47, "", 5, 12, false, 0);
//    test_sites(sites, 48, "", 6, 12, false, 0);
//    test_sites(sites, 49, "", 7, 12, false, 0);
//    test_sites(sites, 50, "", 8, 12, false, 0);
//    test_sites(sites, 51, "", 9, 12, false, 0);
//    test_sites(sites, 52, "", 10, 12, false, 0);
//    test_sites(sites, 53, "", 11, 12, false, 0);
//    test_sites(sites, 54, "", 12, 12, false, 0);
//    test_sites(sites, 55, "", 13, 12, false, 0);
//    test_sites(sites, 56, "", 14, 12, false, 0);
//    test_sites(sites, 57, "", 15, 12, false, 0);
//    test_sites(sites, 58, "", 16, 12, false, 0);
//    test_sites(sites, 59, "", 17, 12, false, 0);
//    test_sites(sites, 60, "", 18, 12, false, 0);
//    test_sites(sites, 61, "", 19, 12, false, 0);
//    test_sites(sites, 62, "", 20, 12, false, 0);
//    test_sites(sites, 63, "", 21, 12, false, 0);
//    test_sites(sites, 64, "", 22, 12, false, 0);
//    test_sites(sites, 65, "", 23, 12, false, 0);
//    test_sites(sites, 66, "", 24, 12, false, 0);
//    test_sites(sites, 67, "", 25, 12, false, 0);
//    test_sites(sites, 68, "", 26, 12, false, 0);
//    test_sites(sites, 69, "", 27, 12, false, 0);
//    test_sites(sites, 70, "", 28, 12, false, 0);
//    test_sites(sites, 71, "", 29, 12, false, 0);
//    test_sites(sites, 72, "", 30, 12, false, 0);
//    test_sites(sites, 73, "", 31, 12, false, 0);
//    test_sites(sites, 74, "", 32, 12, false, 0);
//    test_sites(sites, 75, "", 33, 12, false, 0);
//    test_sites(sites, 76, "", 34, 12, false, 0);
//    test_sites(sites, 77, "", 35, 12, false, 0);
//    test_sites(sites, 78, "", 36, 12, false, 0);
//    test_sites(sites, 79, "", 37, 12, false, 0);
//    test_sites(sites, 80, "", 38, 12, false, 0);
//    test_sites(sites, 81, "", 39, 12, false, 0);
}
}
