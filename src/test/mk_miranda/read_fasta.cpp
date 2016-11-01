#include<string>
#include "gtest/gtest.h"
#include "mr3_core.hpp"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"

namespace {

    class fasta_test : public ::testing::Test
    {
    protected:
        virtual void SetUp() {
            dfile = STRINGIZE(TEST_DATA_PATH);
            seqan::clear(mirna_ids);
            seqan::clear(mrna_ids);
            seqan::clear(mirna_seqs);
            seqan::clear(mrna_seqs);
        }

        virtual void TearDown() {
        }

        int fread_res;
        seqan::CharString dfile;
        mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet mirna_ids;
        mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet mrna_ids;
        mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet mirna_seqs;
        mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet mrna_seqs;
    };

    TEST_F(fasta_test, mirna_fasta) {
        seqan::CharString ifile_mirna;
        ifile_mirna = dfile;
        ifile_mirna += "mir_002.fasta";

        seqan::CharString ifile_mrna;
        ifile_mrna = dfile;
        ifile_mrna += "utr3_001.fasta";

        mr3as::MR3Options options;
        options.mMiRNAFasta = ifile_mirna;
        options.mMRNAFasta = ifile_mrna;

        mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
        coreInput.init_from_args(options);
        int retVal = coreInput.load_seq_from_file();
        EXPECT_EQ(retVal, 0);

        mirna_ids = coreInput.get_mirna_ids();
        EXPECT_EQ(1u, length(mirna_ids));
        EXPECT_STREQ("hsa-miR-0001*MIMAT0000001", seqan::toCString(mirna_ids[0]));

        mirna_seqs = coreInput.get_mirna_seqs();
        EXPECT_EQ(1u, length(mirna_seqs));
        const char *seq1 = "CGUGCCACCCUUUUCCCCAG";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mirna_seqs[0]));
    }

    TEST_F(fasta_test, mrna_fasta) {
        seqan::CharString ifile_mirna;
        ifile_mirna = dfile;
        ifile_mirna += "mir_002.fasta";

        seqan::CharString ifile_mrna;
        ifile_mrna = dfile;
        ifile_mrna += "utr3_011.fasta";

        mr3as::MR3Options options;
        options.mMiRNAFasta = ifile_mirna;
        options.mMRNAFasta = ifile_mrna;

        mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
        coreInput.init_from_args(options);
        int retVal = coreInput.load_seq_from_file();
        EXPECT_EQ(retVal, 0);

        mrna_ids = coreInput.get_mrna_ids();
        EXPECT_EQ(1u, length(mrna_ids));
        EXPECT_STREQ("NM_000001*chr1*-*0000001*9000001", seqan::toCString(mrna_ids[0]));

        mrna_seqs = coreInput.get_mrna_seqs();
        EXPECT_EQ(1u, length(mrna_seqs));
        const char *text = "UUUGAGUUCUGGCCACCAACAAUUUAGUCAUAUCUGAUAGGUACAAAAGA"
                "AAACCAAGAUUUUGAUAUGACCACCUUUCAACACUUUUUACUGCAACUAG";
        EXPECT_STREQ(text, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }

    TEST_F(fasta_test, comp_oput) {
        seqan::CharString dfile1;
        seqan::CharString dfile2;

        dfile1 = dfile;
        dfile1 += "mir_002.fasta";
        dfile2 = dfile;
        dfile2 += "mir_002.fasta";

        int func_res = gtest_compare_two_files(dfile1, dfile2);
        EXPECT_EQ(0u, func_res);
    }

}