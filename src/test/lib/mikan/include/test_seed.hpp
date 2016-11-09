#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mr3_core.hpp"

template <class TSeedSeqs>
class TestSeed : public TestIOMR3AS
{
protected:

    void test_seed(const char *rnastr, int idx, const char *seq_type, bool effective) {
        seqan::RnaString seedseq = rnastr;

        reverseComplement(seedseq);
        comp_two_rnas(mSeedSeqs.get_seed_seq(idx), seedseq);
        EXPECT_STREQ(seq_type, seqan::toCString((seqan::CharString)mSeedSeqs.get_seed_type(idx)));
        EXPECT_EQ(effective, mSeedSeqs.mEffectiveSeeds[idx]);
    }

    TSeedSeqs mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedDef;
};

typedef TestSeed<mr3as::MR3SeedSeqs<seqan::RnaString> > TestSeedMR3AS;
