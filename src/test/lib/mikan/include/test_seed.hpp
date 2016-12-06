#ifndef MIKAN_TEST_SEED_HPP_
#define MIKAN_TEST_SEED_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mr3_core.hpp"

template <class TSeedSeqs, class TTestIO>
class TestSeed : public TTestIO
{
protected:

    void test_seed(const char *rnastr, int idx, const char *seq_type, bool effective, unsigned mmpos) {
        seqan::RnaString seedseq = rnastr;

        reverseComplement(seedseq);
        comp_two_rnas(mSeedSeqs.get_seed_seq(idx), seedseq);
        EXPECT_STREQ(seq_type, seqan::toCString((seqan::CharString)mSeedSeqs.get_seed_type(idx)));
        EXPECT_EQ(effective, mSeedSeqs.mEffectiveSeeds[idx]);
        EXPECT_EQ(mmpos, mSeedSeqs.get_mismatched_pos(idx));
    }

    TSeedSeqs mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedDef;
};

typedef TestSeed<mr3as::MR3SeedSeqs<seqan::RnaString>, TestIOMR3AS> TestSeedMR3AS;

#endif //MIKAN_TEST_SEED_HPP_