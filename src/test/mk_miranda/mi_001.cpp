#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "mr3_core.hpp"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"

namespace {

    class MI001 : public ::testing::Test
    {
    protected:
        virtual void SetUp() {
            dfile = STRINGIZE(TEST_DATA_PATH);
            seqan::CharString ifile_mirna;
            ifile_mirna = dfile;
            ifile_mirna += "mir_001.fasta";

            seqan::CharString ifile_mrna;
            ifile_mrna = dfile;
            ifile_mrna += "utr3_001.fasta";

            options.mMiRNAFasta = ifile_mirna;
            options.mMRNAFasta = ifile_mrna;

            coreInput.init_from_args(options);
            int retVal;
            retVal = coreInput.load_seq_from_file();

            resize(mSeedDef, 6);
            mSeedDef[0] = 'Y';
            mSeedDef[1] = 'Y';
            mSeedDef[2] = 'Y';
            mSeedDef[3] = "+";
            mSeedDef[4] = "1:1";
            mSeedDef[5] = "1";

            seqan::clear(mirna_ids);
            seqan::clear(mrna_ids);
            seqan::clear(mirna_seqs);
            seqan::clear(mrna_seqs);

        }

        virtual void TearDown() {
        }

        int fread_res;
        seqan::CharString dfile;
        mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
        mr3as::MR3Options options;
        mr3as::MR3SeedSeqs<seqan::RnaString> mSeedSeqs;
        seqan::StringSet<seqan::CharString> mSeedDef;
        mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet mirna_ids;
        mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet mrna_ids;
        mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet mirna_seqs;
        mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet mrna_seqs;
    };

    TEST_F(MI001, mirna_id) {
        mirna_ids = coreInput.get_mirna_ids();

        EXPECT_EQ(2u, length(mirna_ids));

        EXPECT_STREQ("hsa-miR-124 MIMAT0000422", seqan::toCString(mirna_ids[0]));
        EXPECT_STREQ("hsa-miR-1 MIMAT0000416", seqan::toCString(mirna_ids[1]));
    }

    TEST_F(MI001, mirna_seq) {
        mirna_seqs = coreInput.get_mirna_seqs();
        EXPECT_EQ(2u, length(mirna_seqs));

        const char *seq1 = "UAAGGCACGCGGUGAAUGCC";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mirna_seqs[0]));

        const char *seq2 = "UGGAAUGUAAAGAAGUAUGUAU";
        EXPECT_STREQ(seq2, seqan::toCString((seqan::CharString)mirna_seqs[1]));
    }

    TEST_F(MI001, get_seed_1) {
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        mSeedDef[1] = 'N';
        mSeedDef[2] = 'N';
        mSeedDef[3] = "0";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "0";

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);

        EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

        seqan::RnaString rnastr = "AAGGCA";
        reverseComplement(rnastr);

        comp_two_rnas(mSeedSeqs.get_seed_seq(0), rnastr);
    }

    TEST_F(MI001, get_seed_6mer) {
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        mSeedDef[1] = 'N';
        mSeedDef[2] = 'N';
        mSeedDef[3] = "0";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "0";

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);

        EXPECT_EQ(1u, length(mSeedSeqs.mEffectiveSeeds));

        seqan::RnaString rnastr = "AAGGCA";
        reverseComplement(rnastr);

        comp_two_rnas(mSeedSeqs.get_seed_seq(0), rnastr);
        EXPECT_TRUE(mSeedSeqs.mEffectiveSeeds[0]);

        const char *seq_type = "6mer";
        EXPECT_STREQ(seq_type, seqan::toCString((seqan::CharString)mSeedSeqs.get_seed_type(0)));
    }

    TEST_F(MI001, get_seed_gut) {
        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        mSeedDef[3] = "1";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "0";

        int n = mSeedSeqs.create_seed_seqs(mSeedDef);
        EXPECT_EQ(0, n);

        EXPECT_EQ(3u, length(mSeedSeqs.mEffectiveSeeds));

        seqan::RnaString rnastr = "AAGGCA";
        reverseComplement(rnastr);
        comp_two_rnas(mSeedSeqs.get_seed_seq(0), rnastr);
        const char *seq_type1 = "6mer";
        EXPECT_STREQ(seq_type1, seqan::toCString((seqan::CharString)mSeedSeqs.get_seed_type(0)));

        rnastr = "AAGACA";
        reverseComplement(rnastr);
        comp_two_rnas(mSeedSeqs.get_seed_seq(1), rnastr);
        const char *seq_type2 = "GUT";
        EXPECT_STREQ(seq_type2, seqan::toCString((seqan::CharString)mSeedSeqs.get_seed_type(1)));

        rnastr = "AAAGCA";
        reverseComplement(rnastr);
        comp_two_rnas(mSeedSeqs.get_seed_seq(2), rnastr);
        const char *seq_type3 = "GUT";
        EXPECT_STREQ(seq_type3, seqan::toCString((seqan::CharString)mSeedSeqs.get_seed_type(2)));
    }
}
