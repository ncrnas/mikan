#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site02GU2 : public TestSiteRH2 {
protected:
    Site02GU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_2.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mGU1";
        mOverlapDef = "orig";
    }

};

TEST_F(Site02GU2, mir1_6mer_gu) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(40u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "6mer", 3, 24, true, 0);

    test_sites(sites, 4, "6mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "6mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "6mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "6mer_GUT", 7, 24, true, 0);

    test_sites(sites, 8, "6mer_GUT", 8, 24, true, 0);
    test_sites(sites, 9, "6mer_GUT", 9, 24, true, 0);
    test_sites(sites, 10, "6mer_GUT", 10, 24, true, 0);
    test_sites(sites, 11, "6mer_GUT", 11, 24, true, 0);
    test_sites(sites, 12, "6mer_GUT", 12, 24, true, 0);

    test_sites(sites, 13, "6mer_GUM", 13, 24, true, 0);
    test_sites(sites, 14, "6mer_GUM", 14, 24, true, 0);
    test_sites(sites, 15, "6mer_GUM", 15, 24, true, 0);
    test_sites(sites, 16, "6mer_GUM", 16, 24, true, 0);

    test_sites(sites, 17, "6mer_GUM", 17, 24, true, 0);
    test_sites(sites, 18, "6mer_GUM", 18, 24, true, 0);
    test_sites(sites, 19, "6mer_GUM", 19, 24, true, 0);
    test_sites(sites, 20, "6mer_GUM", 20, 24, true, 0);
    test_sites(sites, 21, "6mer_GUM", 21, 24, true, 0);

    test_sites(sites, 22, "6mer_GUT", 22, 24, true, 0);
    test_sites(sites, 23, "6mer_GUT", 23, 24, true, 0);
    test_sites(sites, 24, "6mer_GUT", 24, 24, true, 0);
    test_sites(sites, 25, "6mer_GUT", 25, 24, true, 0);

    test_sites(sites, 26, "6mer_GUT", 26, 24, true, 0);
    test_sites(sites, 27, "6mer_GUT", 27, 24, true, 0);
    test_sites(sites, 28, "6mer_GUT", 28, 24, true, 0);
    test_sites(sites, 29, "6mer_GUT", 29, 24, true, 0);
    test_sites(sites, 30, "6mer_GUT", 30, 24, true, 0);

    test_sites(sites, 31, "6mer_GUT", 31, 24, true, 0);
    test_sites(sites, 32, "6mer_GUT", 32, 24, true, 0);
    test_sites(sites, 33, "6mer_GUT", 33, 24, true, 0);
    test_sites(sites, 34, "6mer_GUT", 34, 24, true, 0);

    test_sites(sites, 35, "6mer_GUT", 35, 24, true, 0);
    test_sites(sites, 36, "6mer_GUT", 36, 24, true, 0);
    test_sites(sites, 37, "6mer_GUT", 37, 24, true, 0);
    test_sites(sites, 38, "6mer_GUT", 38, 24, true, 0);
    test_sites(sites, 39, "6mer_GUT", 39, 24, true, 0);
}

TEST_F(Site02GU2, mir1_6mer_gu_plus) {
    mSeedDef[0] = "6mGU+";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(40u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "6mer", 3, 24, true, 0);

    test_sites(sites, 4, "6mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "6mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "6mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "6mer_GUT", 7, 24, true, 0);

    test_sites(sites, 8, "6mer_GUT", 8, 24, true, 0);
    test_sites(sites, 9, "6mer_GUT", 9, 24, true, 0);
    test_sites(sites, 10, "6mer_GUT", 10, 24, true, 0);
    test_sites(sites, 11, "6mer_GUT", 11, 24, true, 0);
    test_sites(sites, 12, "6mer_GUT", 12, 24, true, 0);

    test_sites(sites, 13, "6mer_GUM", 13, 24, true, 0);
    test_sites(sites, 14, "6mer_GUM", 14, 24, true, 0);
    test_sites(sites, 15, "6mer_GUM", 15, 24, true, 0);
    test_sites(sites, 16, "6mer_GUM", 16, 24, true, 0);

    test_sites(sites, 17, "6mer_GUM", 17, 24, true, 0);
    test_sites(sites, 18, "6mer_GUM", 18, 24, true, 0);
    test_sites(sites, 19, "6mer_GUM", 19, 24, true, 0);
    test_sites(sites, 20, "6mer_GUM", 20, 24, true, 0);
    test_sites(sites, 21, "6mer_GUM", 21, 24, true, 0);

    test_sites(sites, 22, "6mer_GUT", 22, 24, true, 0);
    test_sites(sites, 23, "6mer_GUT", 23, 24, true, 0);
    test_sites(sites, 24, "6mer_GUT", 24, 24, true, 0);
    test_sites(sites, 25, "6mer_GUT", 25, 24, true, 0);

    test_sites(sites, 26, "6mer_GUT", 26, 24, true, 0);
    test_sites(sites, 27, "6mer_GUT", 27, 24, true, 0);
    test_sites(sites, 28, "6mer_GUT", 28, 24, true, 0);
    test_sites(sites, 29, "6mer_GUT", 29, 24, true, 0);
    test_sites(sites, 30, "6mer_GUT", 30, 24, true, 0);

    test_sites(sites, 31, "6mer_GUT", 31, 24, true, 0);
    test_sites(sites, 32, "6mer_GUT", 32, 24, true, 0);
    test_sites(sites, 33, "6mer_GUT", 33, 24, true, 0);
    test_sites(sites, 34, "6mer_GUT", 34, 24, true, 0);

    test_sites(sites, 35, "6mer_GUT", 35, 24, true, 0);
    test_sites(sites, 36, "6mer_GUT", 36, 24, true, 0);
    test_sites(sites, 37, "6mer_GUT", 37, 24, true, 0);
    test_sites(sites, 38, "6mer_GUT", 38, 24, true, 0);
    test_sites(sites, 39, "6mer_GUT", 39, 24, true, 0);
}

TEST_F(Site02GU2, mir1_7mer_gu) {
    mSeedDef[0] = "7mGU1";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(20u, sites.get_length());

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, 0);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, 0);
    test_sites(sites, 2, "7mer_GUM", 2, 24, true, 0);
    test_sites(sites, 3, "7mer_GUM", 3, 24, true, 0);
    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "7mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "7mer_GUT", 7, 24, true, 0);
//    test_sites(sites, 8, "7mer_GUT", 8, 24, false);
//    test_sites(sites, 9, "7mer_GUT", 9, 24, false);
//    test_sites(sites, 10, "7mer_GUT", 10, 24, false);
//    test_sites(sites, 11, "7mer_GUT", 11, 24, false);
//    test_sites(sites, 12, "7mer_GUT", 12, 24, false);
    test_sites(sites, 8, "7mer_GUM", 13, 24, true, 0);
    test_sites(sites, 9, "7mer_GUM", 14, 24, true, 0);
    test_sites(sites, 10, "7mer_GUM", 15, 24, true, 0);
    test_sites(sites, 11, "7mer_GUM", 16, 24, true, 0);
//    test_sites(sites, 17, "7mer_GUM", 17, 24, false);
//    test_sites(sites, 18, "7mer_GUM", 18, 24, false);
//    test_sites(sites, 19, "7mer_GUM", 19, 24, false);
//    test_sites(sites, 20, "7mer_GUM", 20, 24, false);
//    test_sites(sites, 21, "7mer_GUM", 21, 24, false);
    test_sites(sites, 12, "7mer_GUT", 22, 24, true, 0);
    test_sites(sites, 13, "7mer_GUT", 23, 24, true, 0);
    test_sites(sites, 14, "7mer_GUT", 24, 24, true, 0);
    test_sites(sites, 15, "7mer_GUT", 25, 24, true, 0);
//    test_sites(sites, 26, "7mer_GUT", 26, 24, false);
//    test_sites(sites, 27, "7mer_GUT", 27, 24, false);
//    test_sites(sites, 28, "7mer_GUT", 28, 24, false);
//    test_sites(sites, 29, "7mer_GUT", 29, 24, false);
//    test_sites(sites, 30, "7mer_GUT", 30, 24, false);
    test_sites(sites, 16, "7mer_GUT", 31, 24, true, 0);
    test_sites(sites, 17, "7mer_GUT", 32, 24, true, 0);
    test_sites(sites, 18, "7mer_GUT", 33, 24, true, 0);
    test_sites(sites, 19, "7mer_GUT", 34, 24, true, 0);
//    test_sites(sites, 35, "7mer_GUT", 35, 24, false);
//    test_sites(sites, 36, "7mer_GUT", 36, 24, false);
//    test_sites(sites, 37, "7mer_GUT", 37, 24, false);
//    test_sites(sites, 38, "7mer_GUT", 38, 24, false);
//    test_sites(sites, 39, "7mer_GUT", 39, 24, false);
}

TEST_F(Site02GU2, mir1_def) {
    mSeedDef[0] = "7mGU+";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(36u, sites.get_length());

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, 0);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, 0);
    test_sites(sites, 2, "7mer_GUM", 2, 24, true, 0);
    test_sites(sites, 3, "7mer_GUM", 3, 24, true, 0);
    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "7mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "7mer_GUT", 7, 24, true, 0);
//    test_sites(sites, 8, "7mer_GUT", 8, 24, false);
    test_sites(sites, 8, "7mer_GUM", 9, 24, true, 0);
    test_sites(sites, 9, "7mer_GUM", 10, 24, true, 0);
    test_sites(sites, 10, "7mer_GUM", 11, 24, true, 0);
    test_sites(sites, 11, "7mer_GUM", 12, 24, true, 0);
    test_sites(sites, 12, "7mer_GUM", 13, 24, true, 0);
    test_sites(sites, 13, "7mer_GUM", 14, 24, true, 0);
    test_sites(sites, 14, "7mer_GUM", 15, 24, true, 0);
    test_sites(sites, 15, "7mer_GUM", 16, 24, true, 0);
//    test_sites(sites, 17, "7mer_GUM", 17, 24, false);
    test_sites(sites, 16, "7mer_GUM", 18, 24, true, 0);
    test_sites(sites, 17, "7mer_GUM", 19, 24, true, 0);
    test_sites(sites, 18, "7mer_GUM", 20, 24, true, 0);
    test_sites(sites, 19, "7mer_GUM", 21, 24, true, 0);
    test_sites(sites, 20, "7mer_GUT", 22, 24, true, 0);
    test_sites(sites, 21, "7mer_GUT", 23, 24, true, 0);
    test_sites(sites, 22, "7mer_GUT", 24, 24, true, 0);
    test_sites(sites, 23, "7mer_GUT", 25, 24, true, 0);
//    test_sites(sites, 26, "7mer_GUT", 26, 24, false);
    test_sites(sites, 24, "7mer_GUM", 27, 24, true, 0);
    test_sites(sites, 25, "7mer_GUM", 28, 24, true, 0);
    test_sites(sites, 26, "7mer_GUM", 29, 24, true, 0);
    test_sites(sites, 27, "7mer_GUM", 30, 24, true, 0);
    test_sites(sites, 28, "7mer_GUT", 31, 24, true, 0);
    test_sites(sites, 29, "7mer_GUT", 32, 24, true, 0);
    test_sites(sites, 30, "7mer_GUT", 33, 24, true, 0);
    test_sites(sites, 31, "7mer_GUT", 34, 24, true, 0);
//    test_sites(sites, 35, "7mer_GUT", 35, 24, false);
    test_sites(sites, 32, "7mer_GUM", 36, 24, true, 0);
    test_sites(sites, 33, "7mer_GUM", 37, 24, true, 0);
    test_sites(sites, 34, "7mer_GUM", 38, 24, true, 0);
    test_sites(sites, 35, "7mer_GUM", 39, 24, true, 0);
}
}
