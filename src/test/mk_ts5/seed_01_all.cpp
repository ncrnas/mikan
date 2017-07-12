#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class SeedAll : public TestSeedTS5 {
protected:
    SeedAll() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_ts5/";
    }
};

TEST_F(SeedAll, mir124_def) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedSeqs.set_mirna_seq(mirna_seqs[0]);
    mikan::TCharSet mNullSet;
    mSeedSeqs.set_flags(mNullSet);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    test_seed3("AAGGCA");
}

TEST_F(SeedAll, mir1_def) {
    read_files();

    mirna_seqs = coreInput.get_mirna_seqs();
    mSeedSeqs.set_mirna_seq(mirna_seqs[1]);
    mikan::TCharSet mNullSet;
    mSeedSeqs.set_flags(mNullSet);
    int n = mSeedSeqs.create_seed_seqs();
    EXPECT_EQ(0, n);
    test_seed3("GGAAUG");
}
}
