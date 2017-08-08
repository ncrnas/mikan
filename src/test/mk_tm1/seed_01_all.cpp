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
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tm1/";
    }
};

TEST_F(SeedAll, mir124_def) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();

    mikan::TCharSet mNullSet;
    mSeedSeqs.set_seed_type_def(mNullSet);
    int n = mSeedSeqs.create_seed_seqs(mirna_seqs[0]);
    EXPECT_EQ(0, n);
    EXPECT_EQ(12u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("AAGGCA", 0, "6mer", true);
    test_seed2("AAGACA", 1, "GUT", true);
    test_seed2("AAAGCA", 2, "GUT", true);

    test_seed2("UAGGCA", 3, "MM", true);
    test_seed2("GAGGCA", 4, "MM", true);
    test_seed2("CAGGCA", 5, "MM", true);

    test_seed2("UAGACA", 6, "GUMM", true);
    test_seed2("GAGACA", 7, "GUMM", true);
    test_seed2("CAGACA", 8, "GUMM", true);
    test_seed2("UAAGCA", 9, "GUMM", true);
    test_seed2("GAAGCA", 10, "GUMM", true);
    test_seed2("CAAGCA", 11, "GUMM", true);
}

TEST_F(SeedAll, mir1_def) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();

    mikan::TCharSet mNullSet;
    mSeedSeqs.set_seed_type_def(mNullSet);
    int n = mSeedSeqs.create_seed_seqs(mirna_seqs[1]);
    EXPECT_EQ(0, n);
    EXPECT_EQ(19u, length(mSeedSeqs.mEffectiveSeeds));

    test_seed2("GGAAUG", 0, "6mer", true);
    test_seed2("GGAAUA", 1, "GUT", true);
    test_seed2("GGAACG", 2, "GUM", true);
    test_seed2("GAAAUG", 3, "GUT", true);
    test_seed2("AGAAUG", 4, "GUT", true);

    test_seed2("UGAAUG", 5, "MM", false);
    test_seed2("CGAAUG", 6, "MM", false);

    test_seed2("UGAAUA", 7, "GUMM", true);
    test_seed2("CGAAUA", 8, "GUMM", true);
    test_seed2("AGAAUA", 9, "GUGU", true);

    test_seed2("UGAACG", 10, "GUMM", true);
    test_seed2("CGAACG", 11, "GUMM", true);
    test_seed2("AGAACG", 12, "GUGU", true);

    test_seed2("UAAAUG", 13, "GUMM", true);
    test_seed2("CAAAUG", 14, "GUMM", true);
    test_seed2("AAAAUG", 15, "GUGU", true);

    test_seed2("UGAAUG", 16, "GUMM", true);
    test_seed2("GGAAUG", 17, "GUMM", false);
    test_seed2("CGAAUG", 18, "GUMM", true);
}

}
