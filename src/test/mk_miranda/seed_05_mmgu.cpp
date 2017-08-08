#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class SeedMMGU : public TestSeedMR3AS {
protected:
    SeedMMGU() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_miranda/";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "1";
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "0";
    }
};

TEST_F(SeedMMGU, mir124_mmgu) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();

    mSeedSeqs.set_seed_type_def(mSeedDef);
    mSeedSeqs.set_flags();
    int n = mSeedSeqs.create_seed_seqs(mirna_seqs[0]);
    EXPECT_EQ(0, n);
    EXPECT_EQ(47u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("AAGGCA", 0, "6mer", true, 0);

    test_seed("AAGACA", 1, "GUT", true, 2);
    test_seed("AAAGCA", 2, "GUT", true, 3);

    test_seed("AAGGCG", 3, "MM", true, 0);
    test_seed("AAGGCC", 4, "MM", true, 0);
    test_seed("AAGGCU", 5, "MM", true, 0);
    test_seed("AAGGUA", 6, "MM", true, 1);
    test_seed("AAGGAA", 7, "MM", true, 1);
    test_seed("AAGGGA", 8, "MM", true, 1);
    test_seed("AAGUCA", 9, "MM", true, 2);
    test_seed("AAGCCA", 10, "MM", true, 2);
    test_seed("AAUGCA", 11, "MM", true, 3);
    test_seed("AACGCA", 12, "MM", true, 3);
    test_seed("AGGGCA", 13, "MM", true, 4);
    test_seed("ACGGCA", 14, "MM", true, 4);
    test_seed("AUGGCA", 15, "MM", true, 4);
    test_seed("GAGGCA", 16, "MM", true, 5);
    test_seed("CAGGCA", 17, "MM", true, 5);
    test_seed("UAGGCA", 18, "MM", true, 5);

    test_seed("AAGACG", 19, "MMGU", true, 0);
    test_seed("AAGACC", 20, "MMGU", true, 0);
    test_seed("AAGACU", 21, "MMGU", true, 0);
    test_seed("AAGAUA", 22, "MMGU", true, 1);
    test_seed("AAGAAA", 23, "MMGU", true, 1);
    test_seed("AAGAGA", 24, "MMGU", true, 1);
    test_seed("AAUACA", 25, "MMGU", true, 3);
    test_seed("AACACA", 26, "MMGU", true, 3);
    test_seed("AGGACA", 27, "MMGU", true, 4);
    test_seed("ACGACA", 28, "MMGU", true, 4);
    test_seed("AUGACA", 29, "MMGU", true, 4);
    test_seed("GAGACA", 30, "MMGU", true, 5);
    test_seed("CAGACA", 31, "MMGU", true, 5);
    test_seed("UAGACA", 32, "MMGU", true, 5);

    test_seed("AAAGCG", 33, "MMGU", true, 0);
    test_seed("AAAGCC", 34, "MMGU", true, 0);
    test_seed("AAAGCU", 35, "MMGU", true, 0);
    test_seed("AAAGUA", 36, "MMGU", true, 1);
    test_seed("AAAGAA", 37, "MMGU", true, 1);
    test_seed("AAAGGA", 38, "MMGU", true, 1);
    test_seed("AAAUCA", 39, "MMGU", true, 2);
    test_seed("AAACCA", 40, "MMGU", true, 2);
    test_seed("AGAGCA", 41, "MMGU", true, 4);
    test_seed("ACAGCA", 42, "MMGU", true, 4);
    test_seed("AUAGCA", 43, "MMGU", true, 4);
    test_seed("GAAGCA", 44, "MMGU", true, 5);
    test_seed("CAAGCA", 45, "MMGU", true, 5);
    test_seed("UAAGCA", 46, "MMGU", true, 5);
}

TEST_F(SeedMMGU, mir1_mmgu) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();

    mSeedSeqs.set_seed_type_def(mSeedDef);
    mSeedSeqs.set_flags();
    int n = mSeedSeqs.create_seed_seqs(mirna_seqs[1]);
    EXPECT_EQ(0, n);
    EXPECT_EQ(67u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("GGAAUG", 0, "6mer", true, 0);

    test_seed("GGAAUA", 1, "GUT", true, 0);
    test_seed("GGAACG", 2, "GUM", true, 1);
    test_seed("GAAAUG", 3, "GUT", true, 4);
    test_seed("AGAAUG", 4, "GUT", true, 5);

    test_seed("GGAAUU", 5, "MM", true, 0);
    test_seed("GGAAUC", 6, "MM", true, 0);
    test_seed("GGAAGG", 7, "MM", true, 1);
    test_seed("GGAAAG", 8, "MM", true, 1);
    test_seed("GGAGUG", 9, "MM", true, 2);
    test_seed("GGACUG", 10, "MM", true, 2);
    test_seed("GGAUUG", 11, "MM", true, 2);
    test_seed("GGGAUG", 12, "MM", true, 3);
    test_seed("GGCAUG", 13, "MM", true, 3);
    test_seed("GGUAUG", 14, "MM", true, 3);
    test_seed("GUAAUG", 15, "MM", true, 4);
    test_seed("GCAAUG", 16, "MM", true, 4);
    test_seed("UGAAUG", 17, "MM", true, 5);
    test_seed("CGAAUG", 18, "MM", true, 5);

    test_seed("GGAAGA", 19, "MMGU", true, 1);
    test_seed("GGAAAA", 20, "MMGU", true, 1);
    test_seed("GGAGUA", 21, "MMGU", true, 2);
    test_seed("GGACUA", 22, "MMGU", true, 2);
    test_seed("GGAUUA", 23, "MMGU", true, 2);
    test_seed("GGGAUA", 24, "MMGU", true, 3);
    test_seed("GGCAUA", 25, "MMGU", true, 3);
    test_seed("GGUAUA", 26, "MMGU", true, 3);
    test_seed("GUAAUA", 27, "MMGU", true, 4);
    test_seed("GCAAUA", 28, "MMGU", true, 4);
    test_seed("UGAAUA", 29, "MMGU", true, 5);
    test_seed("CGAAUA", 30, "MMGU", true, 5);

    test_seed("GGAACU", 31, "MMGU", true, 0);
    test_seed("GGAACC", 32, "MMGU", true, 0);
    test_seed("GGAGCG", 33, "MMGU", true, 2);
    test_seed("GGACCG", 34, "MMGU", true, 2);
    test_seed("GGAUCG", 35, "MMGU", true, 2);
    test_seed("GGGACG", 36, "MMGU", true, 3);
    test_seed("GGCACG", 37, "MMGU", true, 3);
    test_seed("GGUACG", 38, "MMGU", true, 3);
    test_seed("GUAACG", 39, "MMGU", true, 4);
    test_seed("GCAACG", 40, "MMGU", true, 4);
    test_seed("UGAACG", 41, "MMGU", true, 5);
    test_seed("CGAACG", 42, "MMGU", true, 5);

    test_seed("GAAAUU", 43, "MMGU", true, 0);
    test_seed("GAAAUC", 44, "MMGU", true, 0);
    test_seed("GAAAGG", 45, "MMGU", true, 1);
    test_seed("GAAAAG", 46, "MMGU", true, 1);
    test_seed("GAAGUG", 47, "MMGU", true, 2);
    test_seed("GAACUG", 48, "MMGU", true, 2);
    test_seed("GAAUUG", 49, "MMGU", true, 2);
    test_seed("GAGAUG", 50, "MMGU", true, 3);
    test_seed("GACAUG", 51, "MMGU", true, 3);
    test_seed("GAUAUG", 52, "MMGU", true, 3);
    test_seed("UAAAUG", 53, "MMGU", true, 5);
    test_seed("CAAAUG", 54, "MMGU", true, 5);

    test_seed("AGAAUU", 55, "MMGU", true, 0);
    test_seed("AGAAUC", 56, "MMGU", true, 0);
    test_seed("AGAAGG", 57, "MMGU", true, 1);
    test_seed("AGAAAG", 58, "MMGU", true, 1);
    test_seed("AGAGUG", 59, "MMGU", true, 2);
    test_seed("AGACUG", 60, "MMGU", true, 2);
    test_seed("AGAUUG", 61, "MMGU", true, 2);
    test_seed("AGGAUG", 62, "MMGU", true, 3);
    test_seed("AGCAUG", 63, "MMGU", true, 3);
    test_seed("AGUAUG", 64, "MMGU", true, 3);
    test_seed("AUAAUG", 65, "MMGU", true, 4);
    test_seed("ACAAUG", 66, "MMGU", true, 4);
}
}
