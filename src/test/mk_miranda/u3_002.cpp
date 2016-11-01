#include <iostream>
#include "gtest/gtest.h"
#include "test_io.hpp"
#include "mr3_core.hpp"

namespace {

    class U3002 : public TestIO
    {
    protected:
        U3002() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_002.fasta";
            O1FNAME1 = (char *)"test_output1_site_2.txt";
            O1FNAME2 = (char *)"test_output1_mrna_2.txt";
            O2FNAME1 = (char *)"test_output2_site_2.txt";
            O2FNAME2 = (char *)"test_output2_mrna_2.txt";
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

    TEST_F(U3002, mrna_fasta) {
        read_files();
        mrna_ids = coreInput.get_mrna_ids();
        mrna_seqs = coreInput.get_mrna_seqs();

        EXPECT_EQ(1u, length(mrna_ids));

        const char *id1 = "hg18_refGene NM_005498 range=chr19:10544348-10544741 5'pad=0 3'pad=0 "
                "revComp=TRUE strand=- repeatMasking=none";
        EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

        EXPECT_EQ(1u, length(mrna_seqs));
        const char *seq1 = "AAGGGAGAAGAGAUGGGGGCUUGAACACGGGGCUUCCUUACAGCCCCGGA"
                "UGCAGAUUUUAGAGGGAGGGCAGGUGCGGGCUGUGUGUGUCUGUGUGAGG"
                "GCAGGUCCUGGACUUGGCAGUUUCUUGCUCCCAGCACCCGCCCCUUCCUC"
                "ACCUCUUCCUUAUUCCAUAGGCUGGGAGAGAAACUCUCUCUGCUUCCCUC"
                "GCCCUUGGAGCUUUCCCCAUCCCCCUGAUUUUAUAUGAAGAAAUAGAAGA"
                "GGGGCUUGAAGUCCUCCUCGCGAGUGCCUUCUUGCAAUUACCUGCCUUAG"
                "CGGGUGUUGCGGGUCCCUCCUUCACAGCCGCUGAGCCCAGAGGUCCCGCU"
                "GGCCCCUCCUCUGAAUUUUAGGAUGUCAUUAAAAAGAUGAAUCU";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }

    TEST_F(U3002, comp_site) {
        run_main();

        gtest_compare_two_files(o1file1, o2file1);
    }

    TEST_F(U3002, comp_mrna) {
        run_main();

        gtest_compare_two_files(o1file2, o2file2);
    }
}