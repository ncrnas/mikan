#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_main_io.hpp"
#include "test_seed.hpp"

namespace {

    class SDM001 : public TestSeedMR3AS
    {
    protected:
        SDM001() {
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

    TEST_F(SDM001, get_seed_1) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        mSeedDef[1] = 'N';
        mSeedDef[2] = 'N';

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer");
    }

    TEST_F(SDM001, get_seed_2_1) {
        read_files(false);
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedDef[1] = 'N';
        mSeedDef[2] = 'N';

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer");
    }

    TEST_F(SDM001, get_seed_2_2) {
        read_files(false);
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedDef[1] = 'N';
        mSeedDef[2] = 'N';

        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("GGAAUG", 0, "6mer");
    }

    TEST_F(SDM001, get_seed_6mer) {
        read_files(false);
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedDef[1] = 'N';
        mSeedDef[2] = 'N';

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer");
    }

    TEST_F(SDM001, get_seed_gut) {
        read_files(false);
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedDef[3] = "1";

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(3u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer");
        test_seed("AAGACA", 1, "GUT");
        test_seed("AAAGCA", 2, "GUT");
    }

    TEST_F(SDM001, get_seed_gum) {
        read_files(false);
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedDef[3] = "1";

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(3u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer");
        test_seed("AAGACA", 1, "GUT");
        test_seed("AAAGCA", 2, "GUT");
    }
}
