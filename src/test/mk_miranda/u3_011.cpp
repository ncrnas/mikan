#include <iostream>
#include "gtest/gtest.h"
#include "test_io.hpp"
#include "mr3_core.hpp"

namespace {

    class U3011 : public TestIO
    {
    protected:
        U3011() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_011.fasta";
            OFNAME1 = (char *)"test_output_1.txt";
            OFNAME2 = (char *)"test_ds1.txt";
            OMPATH = (char *)"mkmiranda/";
        }

        void read_files() {
            options.mMiRNAFasta = ifile1;
            options.mMRNAFasta = ifile2;

            coreInput.init_from_args(options);
            int retVal = coreInput.load_seq_from_file();
            EXPECT_EQ(retVal, 0);
        }
        mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
        mr3as::MR3Options options;
    };

    TEST_F(U3011, mrna_fasta) {
        read_files();
        mrna_ids = coreInput.get_mrna_ids();
        mrna_seqs = coreInput.get_mrna_seqs();

        EXPECT_EQ(1u, length(mrna_ids));

        const char *id1 = "NM_000001*chr1*-*0000001*9000001";
        EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

        EXPECT_EQ(1u, length(mrna_seqs));
        const char *seq1 = "UUUGAGUUCUGGCCACCAACAAUUUAGUCAUAUCUGAUAGGUACAAAAGA"
                "AAACCAAGAUUUUGAUAUGACCACCUUUCAACACUUUUUACUGCAACUAG";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }
}