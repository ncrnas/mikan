#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

    class SeedNmer : public TestSeedTM1
    {
    protected:
        SeedNmer() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_001.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_tm1/";
        }
    };

    TEST_F(SeedNmer, mir124) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs();
        EXPECT_EQ(0, n);
        EXPECT_EQ(12u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed2("AAGGCA", 0, "6mer", true);
        test_seed2("GAGGCA", 1, "MM", true);
        test_seed2("CAGGCA", 2, "MM", true);
        test_seed2("UAGGCA", 3, "MM", true);

        test_seed2("AAGACA", 4, "GUT", true);
        test_seed2("GAGACA", 5, "GUMM", true);
        test_seed2("CAGACA", 6, "GUMM", true);
        test_seed2("UAGACA", 7, "GUMM", true);

        test_seed2("AAAGCA", 8, "GUT", true);
        test_seed2("GAAGCA", 9, "GUMM", true);
        test_seed2("CAAGCA", 10, "GUMM", true);
        test_seed2("UAAGCA", 11, "GUMM", true);
    }

    TEST_F(SeedNmer, mir1) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs();
        EXPECT_EQ(0, n);
        EXPECT_EQ(20u, length(mSeedSeqs.mEffectiveSeeds));

        test_seed2("GGAAUG", 0, "6mer", true);
        test_seed2("UGAAUG", 1, "MM", true);
        test_seed2("AGAAUG", 2, "GUT", true);
        test_seed2("CGAAUG", 3, "MM", true);

        test_seed2("GGAAUA", 4, "GUT", true);
        test_seed2("UGAAUA", 5, "GUMM", true);
        test_seed2("AGAAUA", 6, "GUGU", true);
        test_seed2("CGAAUA", 7, "GUMM", true);

        test_seed2("GGAACG", 8, "GUM", true);
        test_seed2("UGAACG", 9, "GUMM", true);
        test_seed2("AGAACG", 10, "GUGU", true);
        test_seed2("CGAACG", 11, "GUMM", true);

        test_seed2("GAAAUG", 12, "GUT", true);
        test_seed2("UAAAUG", 13, "GUMM", true);
        test_seed2("AAAAUG", 14, "GUGU", true);
        test_seed2("CAAAUG", 15, "GUMM", true);

        test_seed2("AGAAUG", 16, "GUT", false);
        test_seed2("GGAAUG", 17, "GUMM", false);
        test_seed2("CGAAUG", 18, "GUMM", false);
        test_seed2("UGAAUG", 19, "GUMM", false);
    }
}
