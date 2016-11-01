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
            O1FNAME1 = (char *)"test_output1_site_11.txt";
            O1FNAME2 = (char *)"test_output1_mrna_11.txt";
            O2FNAME1 = (char *)"test_output2_site_11.txt";
            O2FNAME2 = (char *)"test_output2_mrna_11.txt";
            OMPATH = (char *)"mk_miranda/";
        }

        void read_files() {
            (void)options.parseCommandLine(argc, (const char **)argv);
            coreInput.init_from_args(options);
            (void)coreInput.load_seq_from_file();
        }

        void run_main() {
            (void)mr3as::MR3CoreMain(argc, (const char **)argv);
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


    TEST_F(U3011, comp_site) {
        run_main();

        gtest_compare_two_files(o1file1, o2file1);
    }

    TEST_F(U3011, comp_mrna) {
        run_main();

        gtest_compare_two_files(o1file2, o2file2);
    }
}