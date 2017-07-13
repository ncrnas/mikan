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
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        resize(mSeedDef, 1);
        mSeedDef[0] = "7mGU+";
        mOverlapDef = "orig";
    }
};

TEST_F(SeedAll, mir124_def) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedSeqs.set_mirna_seq(mirna_seqs[0]);
    mSeedSeqs.set_flags(mSeedDef);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    EXPECT_EQ(4u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("AAGGCA", 0, "7mer", true);

    test_seed2("AAGACA", 1, "7mer_GUT", true);

    test_seed2("AAAGCA", 2, "7mer_GUT", true);
    test_seed2("AAAACA", 3, "7mer_GUT", true);
}

TEST_F(SeedAll, mir1_def) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedSeqs.set_mirna_seq(mirna_seqs[1]);
    mSeedSeqs.set_flags(mSeedDef);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    EXPECT_EQ(16u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("GGAAUG", 0, "7mer", true);

    test_seed2("GGAAUA", 1, "7mer_GUT", true);

    test_seed2("GGAACG", 2, "7mer_GUM", true);
    test_seed2("GAAAUG", 3, "7mer_GUT", true);

    test_seed2("AGAAUG", 4, "7mer_GUT", true);
    test_seed2("GGAACA", 5, "7mer_GUT", true);
    test_seed2("GAAAUA", 6, "7mer_GUT", true);
    test_seed2("AGAAUA", 7, "7mer_GUT", true);

    test_seed2("GAAACG", 8, "7mer_GUM", true);
    test_seed2("AGAACG", 9, "7mer_GUM", true);
    test_seed2("GAAACA", 10, "7mer_GUM", true);
    test_seed2("AGAACA", 11, "7mer_GUM", true);
    test_seed2("AAAAUG", 12, "7mer_GUT", true);
    test_seed2("AAAAUA", 13, "7mer_GUT", true);
    test_seed2("AAAACG", 14, "7mer_GUT", true);
    test_seed2("AAAACA", 15, "7mer_GUT", true);
}
}
