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
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mer";
        mOverlapDef = "orig";
    }
};

TEST_F(SeedNmer, mir124_6mer) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedSeqs.set_mirna_seq(mirna_seqs[0]);
    mSeedSeqs.set_flags(mSeedDef);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("AAGGCA", 0, "6mer", true);
}

TEST_F(SeedNmer, mir1_6mer) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedSeqs.set_mirna_seq(mirna_seqs[1]);
    mSeedSeqs.set_flags(mSeedDef);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("GGAAUG", 0, "6mer", true);
}

TEST_F(SeedNmer, mir124_7mer) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedDef[0] = "7mer";
    mSeedSeqs.set_mirna_seq(mirna_seqs[0]);
    mSeedSeqs.set_flags(mSeedDef);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("AAGGCA", 0, "7mer", true);
}

TEST_F(SeedNmer, mir1_7mer) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedDef[0] = "7mer";
    mSeedSeqs.set_mirna_seq(mirna_seqs[1]);
    mSeedSeqs.set_flags(mSeedDef);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("GGAAUG", 0, "7mer", true);
}
}
