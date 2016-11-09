#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_fasta.hpp"

namespace {

    class MI001 : public TestFasta
    {
    protected:
        MI001() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_001.fasta";
        }
    };

    TEST_F(MI001, mirna_id) {
        read_files();

        mirna_ids = coreInput.get_mirna_ids();
        EXPECT_EQ(2u, length(mirna_ids));

        EXPECT_STREQ("hsa-miR-124 MIMAT0000422", seqan::toCString(mirna_ids[0]));
        EXPECT_STREQ("hsa-miR-1 MIMAT0000416", seqan::toCString(mirna_ids[1]));
    }

    TEST_F(MI001, mirna_seq) {
        read_files();

        mirna_seqs = coreInput.get_mirna_seqs();
        EXPECT_EQ(2u, length(mirna_seqs));

        const char *seq1 = "UAAGGCACGCGGUGAAUGCC";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mirna_seqs[0]));

        const char *seq2 = "UGGAAUGUAAAGAAGUAUGUAU";
        EXPECT_STREQ(seq2, seqan::toCString((seqan::CharString)mirna_seqs[1]));
    }

}
