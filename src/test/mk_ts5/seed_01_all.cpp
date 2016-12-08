#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

    class SeedNmer : public TestSeedTS5
    {
    protected:
        SeedNmer() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_001.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_ts5/";
        }
    };

    TEST_F(SeedNmer, mir124) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        int n = mSeedSeqs.create_seed_seqs();
        EXPECT_EQ(0, n);
        test_seed3("AAGGCA");
    }

    TEST_F(SeedNmer, mir1) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();
        mSeedSeqs.set_mirna_seq(mirna_seqs[1]);

        int n = mSeedSeqs.create_seed_seqs();
        EXPECT_EQ(0, n);
        test_seed3("GGAAUG");
    }
}
