#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

    class SeedGU : public TestSeedMR3AS
    {
    protected:
        SeedGU() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_001.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_miranda/";

            resize(mSeedDef, 6);
            mSeedDef[0] = 'Y';
            mSeedDef[1] = 'Y';
            mSeedDef[2] = 'Y';
            mSeedDef[3] = "1";
            mSeedDef[4] = "0:0";
            mSeedDef[5] = "0";
        }
    };

    TEST_F(SeedGU, mir124_gu) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(3u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true, 0);
        test_seed("AAGACA", 1, "GUT", true, 2);
        test_seed("AAAGCA", 2, "GUT", true, 3);
    }

    TEST_F(SeedGU, mir1_gu) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(5u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("GGAAUG", 0, "6mer", true, 0);
        test_seed("GGAAUA", 1, "GUT", true, 0);
        test_seed("GGAACG", 2, "GUM", true, 1);
        test_seed("GAAAUG", 3, "GUT", true, 4);
        test_seed("AGAAUG", 4, "GUT", true, 5);
    }

    TEST_F(SeedGU, mir124_gu_plus) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedDef[3] = "+";
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(4u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true, 0);
        test_seed("AAGACA", 1, "GUT", true, 2);
        test_seed("AAAGCA", 2, "GUT", true, 3);

        test_seed("AAAACA", 3, "GU+", true, 3);
    }

    TEST_F(SeedGU, mir1_gu_plus) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedDef[3] = "+";
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(16u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("GGAAUG", 0, "6mer", true, 0);
        test_seed("GGAAUA", 1, "GUT", true, 0);
        test_seed("GGAACG", 2, "GUM", true, 1);
        test_seed("GAAAUG", 3, "GUT", true, 4);
        test_seed("AGAAUG", 4, "GUT", true, 5);

        test_seed("GGAACA", 5, "GU+", true, 1);
        test_seed("GAAAUA", 6, "GU+", true, 4);
        test_seed("AGAAUA", 7, "GU+", true, 5);

        test_seed("GAAACG", 8, "GU+", true, 4);
        test_seed("AGAACG", 9, "GU+", true, 5);
        test_seed("GAAACA", 10, "GU+", true, 4);
        test_seed("AGAACA", 11, "GU+", true, 5);

        test_seed("AAAAUG", 12, "GU+", true, 5);
        test_seed("AAAAUA", 13, "GU+", true, 5);
        test_seed("AAAACG", 14, "GU+", true, 5);
        test_seed("AAAACA", 15, "GU+", true, 5);
    }
}
