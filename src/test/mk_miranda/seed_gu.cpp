#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_main_io.hpp"
#include "test_seed.hpp"

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
            mSeedDef[3] = "0";
            mSeedDef[4] = "0:0";
            mSeedDef[5] = "0";
        }
    };

    TEST_F(SeedGU, get_seed_gut) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedDef[3] = "1";
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(3u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true);
        test_seed("AAGACA", 1, "GUT", true);
        test_seed("AAAGCA", 2, "GUT", true);
    }

    TEST_F(SeedGU, get_seed_gum) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedDef[3] = "1";
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(5u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("GGAAUG", 0, "6mer", true);
        test_seed("GGAAUA", 1, "GUT", true);
        test_seed("GGAACG", 2, "GUM", true);
        test_seed("GAAAUG", 3, "GUT", true);
        test_seed("AGAAUG", 4, "GUT", true);
    }

    TEST_F(SeedGU, get_seed_gu_plus) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedDef[3] = "+";
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(12u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true);
        test_seed("AAGACA", 1, "GUT", true);
        test_seed("AAAGCA", 2, "GUT", true);

        test_seed("AAGACA", 3, "GU+", false);
        test_seed("AAGACA", 4, "GU+", false);
        test_seed("AAAACA", 5, "GU+", true);

        test_seed("AAAGCA", 6, "GU+", false);
        test_seed("AAAACA", 7, "GU+", false);
        test_seed("AAAGCA", 8, "GU+", false);
        test_seed("AAAACA", 9, "GU+", false);
        test_seed("AAAACA", 10, "GU+", false);
        test_seed("AAAACA", 11, "GU+", false);
    }
}
