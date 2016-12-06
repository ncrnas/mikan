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
        EXPECT_EQ(12u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true, 0);
        test_seed("AAGACA", 1, "GUT", true, 2);
        test_seed("AAAGCA", 2, "GUT", true, 3);

        test_seed("AAGACA", 3, "GU+", false, 0);
        test_seed("AAGACA", 4, "GU+", false, 0);
        test_seed("AAAACA", 5, "GU+", true, 0);

        test_seed("AAAGCA", 6, "GU+", false, 0);
        test_seed("AAAACA", 7, "GU+", false, 0);
        test_seed("AAAGCA", 8, "GU+", false, 0);
        test_seed("AAAACA", 9, "GU+", false, 0);
        test_seed("AAAACA", 10, "GU+", false, 0);
        test_seed("AAAACA", 11, "GU+", false, 0);
    }

    TEST_F(SeedGU, mir1_gu_plus) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedDef[3] = "+";
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(80u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("GGAAUG", 0, "6mer", true, 0);
        test_seed("GGAAUA", 1, "GUT", true, 0);
        test_seed("GGAACG", 2, "GUM", true, 1);
        test_seed("GAAAUG", 3, "GUT", true, 4);
        test_seed("AGAAUG", 4, "GUT", true, 5);

        test_seed("GGAAUA", 5, "GU+", false, 0);
        test_seed("GGAAUA", 6, "GU+", false, 0);
        test_seed("GGAACA", 7, "GU+", true, 0);
        test_seed("GAAAUA", 8, "GU+", true, 0);
        test_seed("AGAAUA", 9, "GU+", true, 0);

        test_seed("GGAACG", 10, "GU+", false, 0);
        test_seed("GGAACA", 11, "GU+", false, 0);
        test_seed("GGAACG", 12, "GU+", false, 0);
        test_seed("GAAACG", 13, "GU+", true, 0);
        test_seed("AGAACG", 14, "GU+", true, 0);
        test_seed("GGAACA", 15, "GU+", false, 0);
        test_seed("GGAACA", 16, "GU+", false, 0);
        test_seed("GGAACA", 17, "GU+", false, 0);
        test_seed("GAAACA", 18, "GU+", true, 0);
        test_seed("AGAACA", 19, "GU+", true, 0);

        test_seed("GAAAUG", 20, "GU+", false, 0);
        test_seed("GAAAUA", 21, "GU+", false, 0);
        test_seed("GAAACG", 22, "GU+", false, 0);
        test_seed("GAAAUG", 23, "GU+", false, 0);
        test_seed("AAAAUG", 24, "GU+", true, 0);
        test_seed("GAAAUA", 25, "GU+", false, 0);
        test_seed("GAAAUA", 26, "GU+", false, 0);
        test_seed("GAAACA", 27, "GU+", false, 0);
        test_seed("GAAAUA", 28, "GU+", false, 0);
        test_seed("AAAAUA", 29, "GU+", true, 0);
        test_seed("GAAACG", 30, "GU+", false, 0);
        test_seed("GAAACA", 31, "GU+", false, 0);
        test_seed("GAAACG", 32, "GU+", false, 0);
        test_seed("GAAACG", 33, "GU+", false, 0);
        test_seed("AAAACG", 34, "GU+", true, 0);
        test_seed("GAAACA", 35, "GU+", false, 0);
        test_seed("GAAACA", 36, "GU+", false, 0);
        test_seed("GAAACA", 37, "GU+", false, 0);
        test_seed("GAAACA", 38, "GU+", false, 0);
        test_seed("AAAACA", 39, "GU+", true, 0);

        test_seed("AGAAUG", 40, "GU+", false, 0);
        test_seed("AGAAUA", 41, "GU+", false, 0);
        test_seed("AGAACG", 42, "GU+", false, 0);
        test_seed("AAAAUG", 43, "GU+", false, 0);
        test_seed("AGAAUG", 44, "GU+", false, 0);
        test_seed("AGAAUA", 45, "GU+", false, 0);
        test_seed("AGAAUA", 46, "GU+", false, 0);
        test_seed("AGAACA", 47, "GU+", false, 0);
        test_seed("AAAAUA", 48, "GU+", false, 0);
        test_seed("AGAAUA", 49, "GU+", false, 0);
        test_seed("AGAACG", 50, "GU+", false, 0);
        test_seed("AGAACA", 51, "GU+", false, 0);
        test_seed("AGAACG", 52, "GU+", false, 0);
        test_seed("AAAACG", 53, "GU+", false, 0);
        test_seed("AGAACG", 54, "GU+", false, 0);
        test_seed("AGAACA", 55, "GU+", false, 0);
        test_seed("AGAACA", 56, "GU+", false, 0);
        test_seed("AGAACA", 57, "GU+", false, 0);
        test_seed("AAAACA", 58, "GU+", false, 0);
        test_seed("AGAACA", 59, "GU+", false, 0);
        test_seed("AAAAUG", 60, "GU+", false, 0);
        test_seed("AAAAUA", 61, "GU+", false, 0);
        test_seed("AAAACG", 62, "GU+", false, 0);
        test_seed("AAAAUG", 63, "GU+", false, 0);
        test_seed("AAAAUG", 64, "GU+", false, 0);
        test_seed("AAAAUA", 65, "GU+", false, 0);
        test_seed("AAAAUA", 66, "GU+", false, 0);
        test_seed("AAAACA", 67, "GU+", false, 0);
        test_seed("AAAAUA", 68, "GU+", false, 0);
        test_seed("AAAAUA", 69, "GU+", false, 0);
        test_seed("AAAACG", 70, "GU+", false, 0);
        test_seed("AAAACA", 71, "GU+", false, 0);
        test_seed("AAAACG", 72, "GU+", false, 0);
        test_seed("AAAACG", 73, "GU+", false, 0);
        test_seed("AAAACG", 74, "GU+", false, 0);
        test_seed("AAAACA", 75, "GU+", false, 0);
        test_seed("AAAACA", 76, "GU+", false, 0);
        test_seed("AAAACA", 77, "GU+", false, 0);
        test_seed("AAAACA", 78, "GU+", false, 0);
        test_seed("AAAACA", 79, "GU+", false, 0);
    }
}
