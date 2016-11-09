#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_fasta.hpp"

namespace {

    class MI002 : public TestFasta
    {
    protected:
        MI002() {
            IFNAME1 = (char *)"mir_002.fasta";
            IFNAME2 = (char *)"utr3_001.fasta";
        }
    };

    TEST_F(MI002, mirna_id) {
        read_files();

        mirna_ids = coreInput.get_mirna_ids();
        EXPECT_EQ(1u, length(mirna_ids));

        EXPECT_STREQ("hsa-miR-0001*MIMAT0000001", seqan::toCString(mirna_ids[0]));
    }

    TEST_F(MI002, mirna_seq) {
        read_files();

        mirna_seqs = coreInput.get_mirna_seqs();
        EXPECT_EQ(1u, length(mirna_seqs));

        const char *seq1 = "CGUGCCACCCUUUUCCCCAG";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mirna_seqs[0]));
    }

}
