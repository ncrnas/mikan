#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

    class SeedAll : public TestSeedMR3AS
    {
    protected:
        SeedAll() {
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
            mSeedDef[3] = "+";
            mSeedDef[4] = "1:1";
            mSeedDef[5] = "1";
        }
    };

    TEST_F(SeedAll, mir124_def) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(68u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed("AAGGCA", 0, "6mer", true, 0);

        test_seed("AAGACA", 1, "GUT", true, 2);
        test_seed("AAAGCA", 2, "GUT", true, 3);

        test_seed("AAAACA", 3, "GU+", true, 3);

        test_seed("AAGGCG", 4, "MM", true, 0);
        test_seed("AAGGCC", 5, "MM", true, 0);
        test_seed("AAGGCU", 6, "MM", true, 0);
        test_seed("AAGGUA", 7, "MM", true, 1);
        test_seed("AAGGAA", 8, "MM", true, 1);
        test_seed("AAGGGA", 9, "MM", true, 1);
        test_seed("AAGUCA", 10, "MM", true, 2);
        test_seed("AAGCCA", 11, "MM", true, 2);
        test_seed("AAUGCA", 12, "MM", true, 3);
        test_seed("AACGCA", 13, "MM", true, 3);
        test_seed("AGGGCA", 14, "MM", true, 4);
        test_seed("ACGGCA", 15, "MM", true, 4);
        test_seed("AUGGCA", 16, "MM", true, 4);
        test_seed("GAGGCA", 17, "MM", true, 5);
        test_seed("CAGGCA", 18, "MM", true, 5);
        test_seed("UAGGCA", 19, "MM", true, 5);

        test_seed("AAGACG", 20, "MMGU", true, 0);
        test_seed("AAGACC", 21, "MMGU", true, 0);
        test_seed("AAGACU", 22, "MMGU", true, 0);
        test_seed("AAGAUA", 23, "MMGU", true, 1);
        test_seed("AAGAAA", 24, "MMGU", true, 1);
        test_seed("AAGAGA", 25, "MMGU", true, 1);
        test_seed("AAUACA", 26, "MMGU", true, 3);
        test_seed("AACACA", 27, "MMGU", true, 3);
        test_seed("AGGACA", 28, "MMGU", true, 4);
        test_seed("ACGACA", 29, "MMGU", true, 4);
        test_seed("AUGACA", 30, "MMGU", true, 4);
        test_seed("GAGACA", 31, "MMGU", true, 5);
        test_seed("CAGACA", 32, "MMGU", true, 5);
        test_seed("UAGACA", 33, "MMGU", true, 5);

        test_seed("AAAGCG", 34, "MMGU", true, 0);
        test_seed("AAAGCC", 35, "MMGU", true, 0);
        test_seed("AAAGCU", 36, "MMGU", true, 0);
        test_seed("AAAGUA", 37, "MMGU", true, 1);
        test_seed("AAAGAA", 38, "MMGU", true, 1);
        test_seed("AAAGGA", 39, "MMGU", true, 1);
        test_seed("AAAUCA", 40, "MMGU", true, 2);
        test_seed("AAACCA", 41, "MMGU", true, 2);
        test_seed("AGAGCA", 42, "MMGU", true, 4);
        test_seed("ACAGCA", 43, "MMGU", true, 4);
        test_seed("AUAGCA", 44, "MMGU", true, 4);
        test_seed("GAAGCA", 45, "MMGU", true, 5);
        test_seed("CAAGCA", 46, "MMGU", true, 5);
        test_seed("UAAGCA", 47, "MMGU", true, 5);

        test_seed("AGGCUA", 48, "BT", true, 1);
        test_seed("AGGCGA", 49, "BT", true, 1);
        test_seed("AGGCCA", 50, "BT", true, 1);
        test_seed("AGGCAA", 51, "BT", true, 1);
        test_seed("AGGUCA", 52, "BT", true, 2);
        test_seed("AGGGCA", 53, "BT", true, 2);
        test_seed("AGGCCA", 54, "BT", false, 2);
        test_seed("AGGACA", 55, "BT", true, 2);
        test_seed("AGUGCA", 56, "BT", true, 3);
        test_seed("AGGGCA", 57, "BT", false, 3);
        test_seed("AGCGCA", 58, "BT", true, 3);
        test_seed("AGAGCA", 59, "BT", true, 3);
        test_seed("AUGGCA", 60, "BT", true, 4);
        test_seed("AGGGCA", 61, "BT", false, 4);
        test_seed("ACGGCA", 62, "BT", true, 4);
        test_seed("AAGGCA", 63, "BT", false, 4);
        test_seed("UAGGCA", 64, "BT", true, 5);
        test_seed("GAGGCA", 65, "BT", true, 5);
        test_seed("CAGGCA", 66, "BT", true, 5);
        test_seed("AAGGCA", 67, "BT", false, 5);
    }

    TEST_F(SeedAll, mir1_def) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);
        EXPECT_EQ(98u, length(mSeedSeqs.mEffectiveSeeds));

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

        test_seed("GGAAUU", 16, "MM", true, 0);
        test_seed("GGAAUC", 17, "MM", true, 0);
        test_seed("GGAAGG", 18, "MM", true, 1);
        test_seed("GGAAAG", 19, "MM", true, 1);
        test_seed("GGAGUG", 20, "MM", true, 2);
        test_seed("GGACUG", 21, "MM", true, 2);
        test_seed("GGAUUG", 22, "MM", true, 2);
        test_seed("GGGAUG", 23, "MM", true, 3);
        test_seed("GGCAUG", 24, "MM", true, 3);
        test_seed("GGUAUG", 25, "MM", true, 3);
        test_seed("GUAAUG", 26, "MM", true, 4);
        test_seed("GCAAUG", 27, "MM", true, 4);
        test_seed("UGAAUG", 28, "MM", true, 5);
        test_seed("CGAAUG", 29, "MM", true, 5);

        test_seed("GGAAGA", 30, "MMGU", true, 1);
        test_seed("GGAAAA", 31, "MMGU", true, 1);
        test_seed("GGAGUA", 32, "MMGU", true, 2);
        test_seed("GGACUA", 33, "MMGU", true, 2);
        test_seed("GGAUUA", 34, "MMGU", true, 2);
        test_seed("GGGAUA", 35, "MMGU", true, 3);
        test_seed("GGCAUA", 36, "MMGU", true, 3);
        test_seed("GGUAUA", 37, "MMGU", true, 3);
        test_seed("GUAAUA", 38, "MMGU", true, 4);
        test_seed("GCAAUA", 39, "MMGU", true, 4);
        test_seed("UGAAUA", 40, "MMGU", true, 5);
        test_seed("CGAAUA", 41, "MMGU", true, 5);

        test_seed("GGAACU", 42, "MMGU", true, 0);
        test_seed("GGAACC", 43, "MMGU", true, 0);
        test_seed("GGAGCG", 44, "MMGU", true, 2);
        test_seed("GGACCG", 45, "MMGU", true, 2);
        test_seed("GGAUCG", 46, "MMGU", true, 2);
        test_seed("GGGACG", 47, "MMGU", true, 3);
        test_seed("GGCACG", 48, "MMGU", true, 3);
        test_seed("GGUACG", 49, "MMGU", true, 3);
        test_seed("GUAACG", 50, "MMGU", true, 4);
        test_seed("GCAACG", 51, "MMGU", true, 4);
        test_seed("UGAACG", 52, "MMGU", true, 5);
        test_seed("CGAACG", 53, "MMGU", true, 5);

        test_seed("GAAAUU", 54, "MMGU", true, 0);
        test_seed("GAAAUC", 55, "MMGU", true, 0);
        test_seed("GAAAGG", 56, "MMGU", true, 1);
        test_seed("GAAAAG", 57, "MMGU", true, 1);
        test_seed("GAAGUG", 58, "MMGU", true, 2);
        test_seed("GAACUG", 59, "MMGU", true, 2);
        test_seed("GAAUUG", 60, "MMGU", true, 2);
        test_seed("GAGAUG", 61, "MMGU", true, 3);
        test_seed("GACAUG", 62, "MMGU", true, 3);
        test_seed("GAUAUG", 63, "MMGU", true, 3);
        test_seed("UAAAUG", 64, "MMGU", true, 5);
        test_seed("CAAAUG", 65, "MMGU", true, 5);

        test_seed("AGAAUU", 66, "MMGU", true, 0);
        test_seed("AGAAUC", 67, "MMGU", true, 0);
        test_seed("AGAAGG", 68, "MMGU", true, 1);
        test_seed("AGAAAG", 69, "MMGU", true, 1);
        test_seed("AGAGUG", 70, "MMGU", true, 2);
        test_seed("AGACUG", 71, "MMGU", true, 2);
        test_seed("AGAUUG", 72, "MMGU", true, 2);
        test_seed("AGGAUG", 73, "MMGU", true, 3);
        test_seed("AGCAUG", 74, "MMGU", true, 3);
        test_seed("AGUAUG", 75, "MMGU", true, 3);
        test_seed("AUAAUG", 76, "MMGU", true, 4);
        test_seed("ACAAUG", 77, "MMGU", true, 4);

        test_seed("GAAUUG", 78, "BT", true, 1);
        test_seed("GAAUGG", 79, "BT", true, 1);
        test_seed("GAAUCG", 80, "BT", true, 1);
        test_seed("GAAUAG", 81, "BT", true, 1);
        test_seed("GAAUUG", 82, "BT", false, 2);
        test_seed("GAAGUG", 83, "BT", true, 2);
        test_seed("GAACUG", 84, "BT", true, 2);
        test_seed("GAAAUG", 85, "BT", false, 2);
        test_seed("GAUAUG", 86, "BT", true, 3);
        test_seed("GAGAUG", 87, "BT", true, 3);
        test_seed("GACAUG", 88, "BT", true, 3);
        test_seed("GAAAUG", 89, "BT", false, 3);
        test_seed("GUAAUG", 90, "BT", true, 4);
        test_seed("GGAAUG", 91, "BT", false, 4);
        test_seed("GCAAUG", 92, "BT", true, 4);
        test_seed("GAAAUG", 93, "BT", false, 4);
        test_seed("UGAAUG", 94, "BT", true, 5);
        test_seed("GGAAUG", 95, "BT", false, 5);
        test_seed("CGAAUG", 96, "BT", true, 5);
        test_seed("AGAAUG", 97, "BT", false, 5);
    }
}
