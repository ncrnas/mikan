#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site01Nmer1 : public TestSiteRH2 {
protected:
    Site01Nmer1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mer";
        mOverlapDef = "orig";
    }

};

TEST_F(Site01Nmer1, mir124_6mer) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(53u, sites.get_length());

//    test_sites(sites, 0, "", 0, 0, false);
    test_sites(sites, 0, "6mer", 1, 1, true, 0);
    test_sites(sites, 1, "6mer", 2, 13, true, 0);
    test_sites(sites, 2, "6mer", 3, 14, true, 0);
    test_sites(sites, 3, "6mer", 4, 15, true, 0);
    test_sites(sites, 4, "6mer", 5, 18, true, 0);
    test_sites(sites, 5, "6mer", 6, 19, true, 0);
    test_sites(sites, 6, "6mer", 7, 20, true, 0);
    test_sites(sites, 7, "6mer", 8, 32, true, 0);
    test_sites(sites, 8, "6mer", 9, 33, true, 0);
    test_sites(sites, 9, "6mer", 10, 34, true, 0);
//    test_sites(sites, 11, "", 11, 0, false);
    test_sites(sites, 10, "6mer", 12, 1, true, 0);
    test_sites(sites, 11, "6mer", 13, 13, true, 0);
    test_sites(sites, 12, "6mer", 14, 14, true, 0);
    test_sites(sites, 13, "6mer", 15, 15, true, 0);
    test_sites(sites, 14, "6mer", 16, 18, true, 0);
    test_sites(sites, 15, "6mer", 17, 19, true, 0);
    test_sites(sites, 16, "6mer", 18, 20, true, 0);
    test_sites(sites, 17, "6mer", 19, 32, true, 0);
    test_sites(sites, 18, "6mer", 20, 33, true, 0);
    test_sites(sites, 19, "6mer", 21, 34, true, 0);
    test_sites(sites, 20, "6mer", 22, 1, true, 0);
    test_sites(sites, 21, "6mer", 23, 13, true, 0);
    test_sites(sites, 22, "6mer", 24, 14, true, 0);
    test_sites(sites, 23, "6mer", 25, 15, true, 0);
    test_sites(sites, 24, "6mer", 26, 18, true, 0);
    test_sites(sites, 25, "6mer", 27, 19, true, 0);
    test_sites(sites, 26, "6mer", 28, 20, true, 0);
    test_sites(sites, 27, "6mer", 29, 21, true, 0);
    test_sites(sites, 28, "6mer", 30, 32, true, 0);
    test_sites(sites, 29, "6mer", 31, 33, true, 0);
    test_sites(sites, 30, "6mer", 32, 34, true, 0);
    test_sites(sites, 31, "6mer", 33, 1, true, 0);
    test_sites(sites, 32, "6mer", 34, 13, true, 0);
    test_sites(sites, 33, "6mer", 35, 14, true, 0);
    test_sites(sites, 34, "6mer", 36, 15, true, 0);
    test_sites(sites, 35, "6mer", 37, 18, true, 0);
    test_sites(sites, 36, "6mer", 38, 19, true, 0);
    test_sites(sites, 37, "6mer", 39, 20, true, 0);
    test_sites(sites, 38, "6mer", 40, 21, true, 0);
    test_sites(sites, 39, "6mer", 41, 32, true, 0);
    test_sites(sites, 40, "6mer", 42, 33, true, 0);
    test_sites(sites, 41, "6mer", 43, 34, true, 0);
    test_sites(sites, 42, "6mer", 44, 2, true, 0);
    test_sites(sites, 43, "6mer", 45, 13, true, 0);
    test_sites(sites, 44, "6mer", 46, 14, true, 0);
    test_sites(sites, 45, "6mer", 47, 15, true, 0);
    test_sites(sites, 46, "6mer", 48, 18, true, 0);
    test_sites(sites, 47, "6mer", 49, 19, true, 0);
    test_sites(sites, 48, "6mer", 50, 20, true, 0);
    test_sites(sites, 49, "6mer", 51, 21, true, 0);
    test_sites(sites, 50, "6mer", 52, 32, true, 0);
    test_sites(sites, 51, "6mer", 53, 33, true, 0);
    test_sites(sites, 52, "6mer", 54, 34, true, 0);
}

TEST_F(Site01Nmer1, mir124_7mer) {
    mSeedDef[0] = "7mer";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(33u, sites.get_length());

//    test_sites(sites, 0, "", 0, 0, false);
//    test_sites(sites, 1, "7mer", 1, 1, false);
//    test_sites(sites, 2, "7mer", 2, 13, false);
//    test_sites(sites, 3, "7mer", 3, 14, false);
//    test_sites(sites, 4, "7mer", 4, 15, false);
//    test_sites(sites, 5, "7mer", 5, 18, false);
//    test_sites(sites, 6, "7mer", 6, 19, false);
//    test_sites(sites, 7, "7mer", 7, 20, false);
//    test_sites(sites, 8, "7mer", 8, 32, false);
//    test_sites(sites, 9, "7mer", 9, 33, false);
//    test_sites(sites, 10, "7mer", 10, 34, false);
//    test_sites(sites, 11, "", 11, 0, false);
//    test_sites(sites, 12, "7mer", 12, 1, false);
//    test_sites(sites, 13, "7mer", 13, 13, false);
//    test_sites(sites, 14, "7mer", 14, 14, false);
//    test_sites(sites, 15, "7mer", 15, 15, false);
//    test_sites(sites, 16, "7mer", 16, 18, false);
//    test_sites(sites, 17, "7mer", 17, 19, false);
//    test_sites(sites, 18, "7mer", 18, 20, false);
//    test_sites(sites, 19, "7mer", 19, 32, false);
//    test_sites(sites, 20, "7mer", 20, 33, false);
//    test_sites(sites, 21, "7mer", 21, 34, false);
    test_sites(sites, 0, "7mer", 22, 1, true, 0);
    test_sites(sites, 1, "7mer", 23, 13, true, 0);
    test_sites(sites, 2, "7mer", 24, 14, true, 0);
    test_sites(sites, 3, "7mer", 25, 15, true, 0);
    test_sites(sites, 4, "7mer", 26, 18, true, 0);
    test_sites(sites, 5, "7mer", 27, 19, true, 0);
    test_sites(sites, 6, "7mer", 28, 20, true, 0);
    test_sites(sites, 7, "7mer", 29, 21, true, 0);
    test_sites(sites, 8, "7mer", 30, 32, true, 0);
    test_sites(sites, 9, "7mer", 31, 33, true, 0);
    test_sites(sites, 10, "7mer", 32, 34, true, 0);
    test_sites(sites, 11, "7mer", 33, 1, true, 0);
    test_sites(sites, 12, "7mer", 34, 13, true, 0);
    test_sites(sites, 13, "7mer", 35, 14, true, 0);
    test_sites(sites, 14, "7mer", 36, 15, true, 0);
    test_sites(sites, 15, "7mer", 37, 18, true, 0);
    test_sites(sites, 16, "7mer", 38, 19, true, 0);
    test_sites(sites, 17, "7mer", 39, 20, true, 0);
    test_sites(sites, 18, "7mer", 40, 21, true, 0);
    test_sites(sites, 19, "7mer", 41, 32, true, 0);
    test_sites(sites, 20, "7mer", 42, 33, true, 0);
    test_sites(sites, 21, "7mer", 43, 34, true, 0);
    test_sites(sites, 22, "7mer", 44, 2, true, 0);
    test_sites(sites, 23, "7mer", 45, 13, true, 0);
    test_sites(sites, 24, "7mer", 46, 14, true, 0);
    test_sites(sites, 25, "7mer", 47, 15, true, 0);
    test_sites(sites, 26, "7mer", 48, 18, true, 0);
    test_sites(sites, 27, "7mer", 49, 19, true, 0);
    test_sites(sites, 28, "7mer", 50, 20, true, 0);
    test_sites(sites, 29, "7mer", 51, 21, true, 0);
    test_sites(sites, 30, "7mer", 52, 32, true, 0);
    test_sites(sites, 31, "7mer", 53, 33, true, 0);
    test_sites(sites, 32, "7mer", 54, 34, true, 0);
}

TEST_F(Site01Nmer1, mir124_def) {
    mSeedDef[0] = "7mGU+";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(33u, sites.get_length());

//    test_sites(sites, 0, "", 0, 0, false);
//    test_sites(sites, 1, "7mer", 1, 1, false);
//    test_sites(sites, 2, "7mer", 2, 13, false);
//    test_sites(sites, 3, "7mer", 3, 14, false);
//    test_sites(sites, 4, "7mer", 4, 15, false);
//    test_sites(sites, 5, "7mer", 5, 18, false);
//    test_sites(sites, 6, "7mer", 6, 19, false);
//    test_sites(sites, 7, "7mer", 7, 20, false);
//    test_sites(sites, 8, "7mer", 8, 32, false);
//    test_sites(sites, 9, "7mer", 9, 33, false);
//    test_sites(sites, 10, "7mer", 10, 34, false);
//    test_sites(sites, 11, "", 11, 0, false);
//    test_sites(sites, 12, "7mer", 12, 1, false);
//    test_sites(sites, 13, "7mer", 13, 13, false);
//    test_sites(sites, 14, "7mer", 14, 14, false);
//    test_sites(sites, 15, "7mer", 15, 15, false);
//    test_sites(sites, 16, "7mer", 16, 18, false);
//    test_sites(sites, 17, "7mer", 17, 19, false);
//    test_sites(sites, 18, "7mer", 18, 20, false);
//    test_sites(sites, 19, "7mer", 19, 32, false);
//    test_sites(sites, 20, "7mer", 20, 33, false);
//    test_sites(sites, 21, "7mer", 21, 34, false);
    test_sites(sites, 0, "7mer", 22, 1, true, 0);
    test_sites(sites, 1, "7mer", 23, 13, true, 0);
    test_sites(sites, 2, "7mer", 24, 14, true, 0);
    test_sites(sites, 3, "7mer", 25, 15, true, 0);
    test_sites(sites, 4, "7mer", 26, 18, true, 0);
    test_sites(sites, 5, "7mer", 27, 19, true, 0);
    test_sites(sites, 6, "7mer", 28, 20, true, 0);
    test_sites(sites, 7, "7mer", 29, 21, true, 0);
    test_sites(sites, 8, "7mer", 30, 32, true, 0);
    test_sites(sites, 9, "7mer", 31, 33, true, 0);
    test_sites(sites, 10, "7mer", 32, 34, true, 0);
    test_sites(sites, 11, "7mer", 33, 1, true, 0);
    test_sites(sites, 12, "7mer", 34, 13, true, 0);
    test_sites(sites, 13, "7mer", 35, 14, true, 0);
    test_sites(sites, 14, "7mer", 36, 15, true, 0);
    test_sites(sites, 15, "7mer", 37, 18, true, 0);
    test_sites(sites, 16, "7mer", 38, 19, true, 0);
    test_sites(sites, 17, "7mer", 39, 20, true, 0);
    test_sites(sites, 18, "7mer", 40, 21, true, 0);
    test_sites(sites, 19, "7mer", 41, 32, true, 0);
    test_sites(sites, 20, "7mer", 42, 33, true, 0);
    test_sites(sites, 21, "7mer", 43, 34, true, 0);
    test_sites(sites, 22, "7mer", 44, 2, true, 0);
    test_sites(sites, 23, "7mer", 45, 13, true, 0);
    test_sites(sites, 24, "7mer", 46, 14, true, 0);
    test_sites(sites, 25, "7mer", 47, 15, true, 0);
    test_sites(sites, 26, "7mer", 48, 18, true, 0);
    test_sites(sites, 27, "7mer", 49, 19, true, 0);
    test_sites(sites, 28, "7mer", 50, 20, true, 0);
    test_sites(sites, 29, "7mer", 51, 21, true, 0);
    test_sites(sites, 30, "7mer", 52, 32, true, 0);
    test_sites(sites, 31, "7mer", 53, 33, true, 0);
    test_sites(sites, 32, "7mer", 54, 34, true, 0);
}

}
