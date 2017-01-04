#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

    class SeedBT : public TestSeedMR3AS
    {
    protected:
        SeedBT() {
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
            mSeedDef[5] = "1";
        }
    };

    TEST_F(SeedBT, mir124_bt) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(21u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true, 0);

        test_seed("AGGCUA", 1, "BT", true, 1);
        test_seed("AGGCGA", 2, "BT", true, 1);
        test_seed("AGGCCA", 3, "BT", true, 1);
        test_seed("AGGCAA", 4, "BT", true, 1);

        test_seed("AGGUCA", 5, "BT", true, 2);
        test_seed("AGGGCA", 6, "BT", true, 2);
        test_seed("AGGCCA", 7, "BT", false, 2);
        test_seed("AGGACA", 8, "BT", true, 2);

        test_seed("AGUGCA", 9, "BT", true, 3);
        test_seed("AGGGCA", 10, "BT", false, 3);
        test_seed("AGCGCA", 11, "BT", true, 3);
        test_seed("AGAGCA", 12, "BT", true, 3);

        test_seed("AUGGCA", 13, "BT", true, 4);
        test_seed("AGGGCA", 14, "BT", false, 4);
        test_seed("ACGGCA", 15, "BT", true, 4);
        test_seed("AAGGCA", 16, "BT", false, 4);

        test_seed("UAGGCA", 17, "BT", true, 5);
        test_seed("GAGGCA", 18, "BT", true, 5);
        test_seed("CAGGCA", 19, "BT", true, 5);
        test_seed("AAGGCA", 20, "BT", false, 5);
    }

    TEST_F(SeedBT, mir1_bt) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(21u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("GGAAUG", 0, "6mer", true, 0);

        test_seed("GAAUUG", 1, "BT", true, 1);
        test_seed("GAAUGG", 2, "BT", true, 1);
        test_seed("GAAUCG", 3, "BT", true, 1);
        test_seed("GAAUAG", 4, "BT", true, 1);

        test_seed("GAAUUG", 5, "BT", false, 2);
        test_seed("GAAGUG", 6, "BT", true, 2);
        test_seed("GAACUG", 7, "BT", true, 2);
        test_seed("GAAAUG", 8, "BT", true, 2);

        test_seed("GAUAUG", 9, "BT", true, 3);
        test_seed("GAGAUG", 10, "BT", true, 3);
        test_seed("GACAUG", 11, "BT", true, 3);
        test_seed("GAAAUG", 12, "BT", false, 3);

        test_seed("GUAAUG", 13, "BT", true, 4);
        test_seed("GGAAUG", 14, "BT", false, 4);
        test_seed("GCAAUG", 15, "BT", true, 4);
        test_seed("GAAAUG", 16, "BT", false, 4);

        test_seed("UGAAUG", 17, "BT", true, 5);
        test_seed("GGAAUG", 18, "BT", false, 5);
        test_seed("CGAAUG", 19, "BT", true, 5);
        test_seed("AGAAUG", 20, "BT", true, 5);
    }
}
