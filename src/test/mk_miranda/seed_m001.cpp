#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_io.hpp"
#include "mr3_core.hpp"

namespace {

    class SDM001 : public TestIOMR3AS
    {
    protected:
        SDM001() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_001.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_miranda/";
        }
    };

    TEST_F(SDM001, get_seed_1) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
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

    TEST_F(SDM001, get_seed_6mer) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
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

    TEST_F(SDM001, get_seed_gut) {
        read_files(false);

        mirna_seqs = coreInput.get_mirna_seqs();

        mSeedSeqs.set_mirna_seq(mirna_seqs[0]);

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
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
