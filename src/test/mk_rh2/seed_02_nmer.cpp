#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class SeedNmer : public TestSeedRH2 {
protected:
    SeedNmer() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mer";
        mOverlapDef = "orig";
    }
};

TEST_F(SeedNmer, mir124_6mer) {
    create_seed_seqs(0);
    
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("AAGGCA", 0, "6mer", true, 0);
}

TEST_F(SeedNmer, mir1_6mer) {
    create_seed_seqs(1);
    
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("GGAAUG", 0, "6mer", true, 0);
}

TEST_F(SeedNmer, mir124_7mer) {
    mSeedDef[0] = "7mer";
    create_seed_seqs(0);
    
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("AAGGCA", 0, "7mer", true, 0);
}

TEST_F(SeedNmer, mir1_7mer) {
    mSeedDef[0] = "7mer";
    create_seed_seqs(1);
    
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("GGAAUG", 0, "7mer", true, 0);
}
}
