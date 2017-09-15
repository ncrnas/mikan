#include <iostream>
#include "gtest/gtest.h"
#include "test_fasta.hpp"

namespace {

class U3101 : public TestFasta {
protected:
    U3101() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_101.fasta";
    }
};

TEST_F(U3101, mrna_fasta) {
    read_files();
    mrna_ids = coreInput.get_mrna_ids();
    mrna_seqs = coreInput.get_mrna_seqs();

    EXPECT_EQ(3u, length(mrna_ids));
    EXPECT_EQ(3u, length(mrna_seqs));

    const char *id1 = "RNA5Test1";
    EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

    const char *seq1 = "ANCNNNGNNNNNNNNNNNNUUNNNNNANCNNNGNNNNNNNNNNNNUUNNNNN";
    EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString) mrna_seqs[0]));

    const char *id2 = "RNA5Test2";
    EXPECT_STREQ(id2, seqan::toCString(mrna_ids[1]));

    const char *seq2 = "NNNNNNNNNN";
    EXPECT_STREQ(seq2, seqan::toCString((seqan::CharString) mrna_seqs[1]));

    const char *id3 = "NM_001256435*chr17_random*-*2420360*2422298*part1";
    EXPECT_STREQ(id3, seqan::toCString(mrna_ids[2]));

    const char *seq3 = "UGUCCAGCACAGUGGAGAGGNCCGGAAGUGAGACGGGCAGACGGCACCUG"
            "CAGCCUGAAACGCACCGCUCNUNGCGUGCGCCCCCACCUGGUCCCCGGAU"
            "GCCCCCACCACCUGGACAGANGNCACACUGACUGCCCACCCAGCUGUGGC";
    EXPECT_STREQ(seq3, seqan::toCString((seqan::CharString) mrna_seqs[2]));
}
}
