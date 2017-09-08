#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class SeedAll : public TestSeedRH2 {
protected:
    SeedAll() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "7mGU+";
        mOverlapDef = "orig";
    }
};

TEST_F(SeedAll, mir124_def) {
    create_seed_seqs(0);
    
    EXPECT_EQ(4u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("AAGGCA", 0, "7mer", true, 0);

    test_seed("AAGACA", 1, "7mer_GUT", true, 2);

    test_seed("AAAGCA", 2, "7mer_GUT", true, 3);
    test_seed("AAAACA", 3, "7mer_GUT", true, 3);
}

TEST_F(SeedAll, mir1_def) {
    create_seed_seqs(1);
    
    EXPECT_EQ(16u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("GGAAUG", 0, "7mer", true, 0);

    test_seed("GGAAUA", 1, "7mer_GUT", true, 0);

    test_seed("GGAACG", 2, "7mer_GUM", true, 1);
    test_seed("GAAAUG", 3, "7mer_GUT", true, 4);

    test_seed("AGAAUG", 4, "7mer_GUT", true, 5);
    test_seed("GGAACA", 5, "7mer_GUT", true, 1);
    test_seed("GAAAUA", 6, "7mer_GUT", true, 4);
    test_seed("AGAAUA", 7, "7mer_GUT", true, 5);

    test_seed("GAAACG", 8, "7mer_GUM", true, 4);
    test_seed("AGAACG", 9, "7mer_GUM", true, 5);
    test_seed("GAAACA", 10, "7mer_GUM", true, 4);
    test_seed("AGAACA", 11, "7mer_GUM", true, 5);
    test_seed("AAAAUG", 12, "7mer_GUT", true, 5);
    test_seed("AAAAUA", 13, "7mer_GUT", true, 5);
    test_seed("AAAACG", 14, "7mer_GUT", true, 5);
    test_seed("AAAACA", 15, "7mer_GUT", true, 5);
}
}
