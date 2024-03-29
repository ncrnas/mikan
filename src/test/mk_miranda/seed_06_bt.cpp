#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class SeedBT : public TestSeedMR3AS {
protected:
    SeedBT() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "0";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "1";
    }
};

TEST_F(SeedBT, mir124_bt) {
    create_seed_seqs(0);

    EXPECT_EQ(25u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("AAGGCA", 0, "6mer", true, 0);

    test_seed("AGGCAU", 1, "BT", true, 0);
    test_seed("AGGCAG", 2, "BT", true, 0);
    test_seed("AGGCAC", 3, "BT", true, 0);
    test_seed("AGGCAA", 4, "BT", true, 0);

    test_seed("AGGCUA", 5, "BT", true, 1);
    test_seed("AGGCGA", 6, "BT", true, 1);
    test_seed("AGGCCA", 7, "BT", true, 1);
    test_seed("AGGCAA", 8, "BT", false, 1);

    test_seed("AGGUCA", 9, "BT", true, 2);
    test_seed("AGGGCA", 10, "BT", true, 2);
    test_seed("AGGCCA", 11, "BT", false, 2);
    test_seed("AGGACA", 12, "BT", true, 2);

    test_seed("AGUGCA", 13, "BT", true, 3);
    test_seed("AGGGCA", 14, "BT", false, 3);
    test_seed("AGCGCA", 15, "BT", true, 3);
    test_seed("AGAGCA", 16, "BT", true, 3);

    test_seed("AUGGCA", 17, "BT", true, 4);
    test_seed("AGGGCA", 18, "BT", false, 4);
    test_seed("ACGGCA", 19, "BT", true, 4);
    test_seed("AAGGCA", 20, "BT", false, 4);

    test_seed("UAGGCA", 21, "BT", true, 5);
    test_seed("GAGGCA", 22, "BT", true, 5);
    test_seed("CAGGCA", 23, "BT", true, 5);
    test_seed("AAGGCA", 24, "BT", false, 5);
}

TEST_F(SeedBT, mir1_bt) {
    create_seed_seqs(1);

    EXPECT_EQ(25u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("GGAAUG", 0, "6mer", true, 0);

    test_seed("GAAUGU", 1, "BT", true, 0);
    test_seed("GAAUGG", 2, "BT", true, 0);
    test_seed("GAAUGC", 3, "BT", true, 0);
    test_seed("GAAUGA", 4, "BT", true, 0);

    test_seed("GAAUUG", 5, "BT", true, 1);
    test_seed("GAAUGG", 6, "BT", false, 1);
    test_seed("GAAUCG", 7, "BT", true, 1);
    test_seed("GAAUAG", 8, "BT", true, 1);

    test_seed("GAAUUG", 9, "BT", false, 2);
    test_seed("GAAGUG", 10, "BT", true, 2);
    test_seed("GAACUG", 11, "BT", true, 2);
    test_seed("GAAAUG", 12, "BT", true, 2);

    test_seed("GAUAUG", 13, "BT", true, 3);
    test_seed("GAGAUG", 14, "BT", true, 3);
    test_seed("GACAUG", 15, "BT", true, 3);
    test_seed("GAAAUG", 16, "BT", false, 3);

    test_seed("GUAAUG", 17, "BT", true, 4);
    test_seed("GGAAUG", 18, "BT", false, 4);
    test_seed("GCAAUG", 19, "BT", true, 4);
    test_seed("GAAAUG", 20, "BT", false, 4);

    test_seed("UGAAUG", 21, "BT", true, 5);
    test_seed("GGAAUG", 22, "BT", false, 5);
    test_seed("CGAAUG", 23, "BT", true, 5);
    test_seed("AGAAUG", 24, "BT", true, 5);
}
}
