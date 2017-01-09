#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

    class SeedAll : public TestSeedTSSVM
    {
    protected:
        SeedAll() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_001.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_tssvm/";
        }
    };

    TEST_F(SeedAll, mir124_def) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs();
        EXPECT_EQ(0, n);
        EXPECT_EQ(49u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true, 0);

        test_seed("AAAGCA", 1, "GUM", true, 3);
        test_seed("AAGACA", 2, "GUM", true, 4);

        test_seed("AGGCAU", 3, "BT", true, 7);
        test_seed("AGGCAG", 4, "BT", true, 7);
        test_seed("AGGCAC", 5, "BT", true, 7);
        test_seed("AGGCAA", 6, "BT", true, 7);

        test_seed("AGGCUA", 7, "BT", true, 6);
        test_seed("AGGCGA", 8, "BT", true, 6);
        test_seed("AGGCCA", 9, "BT", true, 6);
        test_seed("AGGCAA", 10, "BT", true, 6);

        test_seed("AGGUCA", 11, "BT", true, 5);
        test_seed("AGGGCA", 12, "BT", true, 5);
        test_seed("AGGCCA", 13, "BT", true, 5);
        test_seed("AGGACA", 14, "BT", true, 5);

        test_seed("AGUGCA", 15, "BT", true, 4);
        test_seed("AGGGCA", 16, "BT", true, 4);
        test_seed("AGCGCA", 17, "BT", true, 4);
        test_seed("AGAGCA", 18, "BT", true, 4);

        test_seed("AUGGCA", 19, "BT", true, 3);
        test_seed("AGGGCA", 20, "BT", true, 3);
        test_seed("ACGGCA", 21, "BT", true, 3);
        test_seed("AAGGCA", 22, "BT", true, 3);

        test_seed("UAGGCA", 23, "BT", true, 2);
        test_seed("GAGGCA", 24, "BT", true, 2);
        test_seed("CAGGCA", 25, "BT", true, 2);
        test_seed("AAGGCA", 26, "BT", true, 2);

        test_seed("AAGGCC", 27, "BM", true, 6);
        test_seed("AAGGAC", 28, "BM", true, 5);
        test_seed("AAGCAC", 29, "BM", true, 4);
        test_seed("AAGCAC", 30, "BM", true, 3);
        test_seed("AGGCAC", 31, "BM", true, 2);
        test_seed("AGGCAC", 32, "BM", true, 1);

        test_seed("CAGGCA", 33, "LP", true, 1);
        test_seed("GAGGCA", 34, "LP", true, 1);
        test_seed("UAGGCA", 35, "LP", true, 1);

        test_seed("ACGGCA", 36, "LP", true, 2);
        test_seed("AGGGCA", 37, "LP", true, 2);
        test_seed("AUGGCA", 38, "LP", true, 2);

        test_seed("AACGCA", 39, "LP", true, 3);
        test_seed("AAUGCA", 40, "LP", true, 3);

        test_seed("AAGCCA", 41, "LP", true, 4);
        test_seed("AAGUCA", 42, "LP", true, 4);

        test_seed("AAGGAA", 43, "LP", true, 5);
        test_seed("AAGGGA", 44, "LP", true, 5);
        test_seed("AAGGUA", 45, "LP", true, 5);

        test_seed("AAGGCC", 46, "LP", true, 6);
        test_seed("AAGGCG", 47, "LP", true, 6);
        test_seed("AAGGCU", 48, "LP", true, 6);
    }

    TEST_F(SeedAll, mir1_def) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs();
        EXPECT_EQ(0, n);
        EXPECT_EQ(49u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("GGAAUG", 0, "6mer", true, 0);

        test_seed("AGAAUG", 1, "GUM", true, 1);
        test_seed("GAAAUG", 2, "GUM", true, 2);
        test_seed("GGAACG", 3, "GUT", true, 5);
        test_seed("GGAAUA", 4, "GUM", true, 6);

        test_seed("GAAUGU", 5, "BT", true, 7);
        test_seed("GAAUGG", 6, "BT", true, 7);
        test_seed("GAAUGC", 7, "BT", true, 7);
        test_seed("GAAUGA", 8, "BT", true, 7);

        test_seed("GAAUUG", 9, "BT", true, 6);
        test_seed("GAAUGG", 10, "BT", true, 6);
        test_seed("GAAUCG", 11, "BT", true, 6);
        test_seed("GAAUAG", 12, "BT", true, 6);

        test_seed("GAAUUG", 13, "BT", true, 5);
        test_seed("GAAGUG", 14, "BT", true, 5);
        test_seed("GAACUG", 15, "BT", true, 5);
        test_seed("GAAAUG", 16, "BT", true, 5);

        test_seed("GAUAUG", 17, "BT", true, 4);
        test_seed("GAGAUG", 18, "BT", true, 4);
        test_seed("GACAUG", 19, "BT", true, 4);
        test_seed("GAAAUG", 20, "BT", true, 4);

        test_seed("GUAAUG", 21, "BT", true, 3);
        test_seed("GGAAUG", 22, "BT", true, 3);
        test_seed("GCAAUG", 23, "BT", true, 3);
        test_seed("GAAAUG", 24, "BT", true, 3);

        test_seed("UGAAUG", 25, "BT", true, 2);
        test_seed("GGAAUG", 26, "BT", true, 2);
        test_seed("CGAAUG", 27, "BT", true, 2);
        test_seed("AGAAUG", 28, "BT", true, 2);

        test_seed("GGAAUU", 29, "BM", true, 6);
        test_seed("GGAAGU", 30, "BM", true, 5);
        test_seed("GGAUGU", 31, "BM", true, 4);
        test_seed("GGAUGU", 32, "BM", true, 3);
        test_seed("GAAUGU", 33, "BM", true, 2);
        test_seed("GAAUGU", 34, "BM", true, 1);

        test_seed("CGAAUG", 35, "LP", true, 1);
        test_seed("UGAAUG", 36, "LP", true, 1);

        test_seed("GCAAUG", 37, "LP", true, 2);
        test_seed("GUAAUG", 38, "LP", true, 2);

        test_seed("GGCAUG", 39, "LP", true, 3);
        test_seed("GGGAUG", 40, "LP", true, 3);
        test_seed("GGUAUG", 41, "LP", true, 3);

        test_seed("GGACUG", 42, "LP", true, 4);
        test_seed("GGAGUG", 43, "LP", true, 4);
        test_seed("GGAUUG", 44, "LP", true, 4);

        test_seed("GGAAAG", 45, "LP", true, 5);
        test_seed("GGAAGG", 46, "LP", true, 5);

        test_seed("GGAAUC", 47, "LP", true, 6);
        test_seed("GGAAUU", 48, "LP", true, 6);
    }
}
