#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class SeedAll : public TestSeedTM1 {
protected:
    SeedAll() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";
    }
};

TEST_F(SeedAll, mir124_def) {
    create_seed_seqs(0);

    EXPECT_EQ(12u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("AAGGCA", 0, "6mer", true, 0);
    test_seed("AAGACA", 1, "GUT", true, 2);
    test_seed("AAAGCA", 2, "GUT", true, 3);

    test_seed("UAGGCA", 3, "MM", true, 5);
    test_seed("GAGGCA", 4, "MM", true, 5);
    test_seed("CAGGCA", 5, "MM", true, 5);

    test_seed("UAGACA", 6, "GUMM", true, 5);
    test_seed("GAGACA", 7, "GUMM", true, 5);
    test_seed("CAGACA", 8, "GUMM", true, 5);
    test_seed("UAAGCA", 9, "GUMM", true, 5);
    test_seed("GAAGCA", 10, "GUMM", true, 5);
    test_seed("CAAGCA", 11, "GUMM", true, 5);
}

TEST_F(SeedAll, mir1_def) {
    create_seed_seqs(1);

    EXPECT_EQ(19u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed("GGAAUG", 0, "6mer", true, 0);
    test_seed("GGAAUA", 1, "GUT", true, 0);
    test_seed("GGAACG", 2, "GUM", true, 1);
    test_seed("GAAAUG", 3, "GUT", true, 4);
    test_seed("AGAAUG", 4, "GUT", true, 5);

    test_seed("UGAAUG", 5, "MM", false, 5);
    test_seed("CGAAUG", 6, "MM", false, 5);

    test_seed("UGAAUA", 7, "GUMM", true, 5);
    test_seed("CGAAUA", 8, "GUMM", true, 5);
    test_seed("AGAAUA", 9, "GUGU", true, 5);

    test_seed("UGAACG", 10, "GUMM", true, 5);
    test_seed("CGAACG", 11, "GUMM", true, 5);
    test_seed("AGAACG", 12, "GUGU", true, 5);

    test_seed("UAAAUG", 13, "GUMM", true, 5);
    test_seed("CAAAUG", 14, "GUMM", true, 5);
    test_seed("AAAAUG", 15, "GUGU", true, 5);

    test_seed("UGAAUG", 16, "GUMM", true, 5);
    test_seed("GGAAUG", 17, "GUMM", false, 5);
    test_seed("CGAAUG", 18, "GUMM", true, 5);
}

}
