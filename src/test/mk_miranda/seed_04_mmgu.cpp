#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_main_io.hpp"
#include "test_seed.hpp"

namespace {

    class SeedMMGU : public TestSeedMR3AS
    {
    protected:
        SeedMMGU() {
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

    TEST_F(SeedMMGU, mir124_mmgu) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedDef[3] = "1";
        mSeedDef[4] = "1:1";
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
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
}
