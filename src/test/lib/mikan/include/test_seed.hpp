#ifndef MIKAN_TEST_SEED_HPP_
#define MIKAN_TEST_SEED_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mk_option.hpp"
#include "test_main_io.hpp"

template<class TOptions, class TSeedSeqs>
class TestSeed : public TestIOBase {
public:

    TestSeed(): mSeedSeqs(mOpts) {}

protected:

    void test_seed(const char *rnastr, int idx, const char *seq_type, bool effective, unsigned mmpos) {
        seqan::RnaString seedseq = rnastr;

        reverseComplement(seedseq);

        if (idx < length(mSeedSeqs.mEffectiveSeeds)) {
            comp_two_rnas(mSeedSeqs.get_seed_seq(idx), seedseq);
            EXPECT_STREQ(seq_type, seqan::toCString((seqan::CharString) mSeedSeqs.get_seed_type(idx)));
            EXPECT_EQ(effective, mSeedSeqs.mEffectiveSeeds[idx]);
            EXPECT_EQ(mmpos, mSeedSeqs.get_mismatched_pos(idx));
        }
    }

    void create_seed_seqs(unsigned pIdx) {
        read_and_set_seqs();

        mSeedSeqs.set_seed_type_def(mSeedDef);
        mSeedSeqs.set_flags();
        int ret = mSeedSeqs.create_seed_seqs(mirna_seqs[pIdx]);
        EXPECT_EQ(0, ret);
    }

    TOptions mOpts;
    TSeedSeqs mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedDef;
    seqan::CharString mOverlapDef;
};

#endif //MIKAN_TEST_SEED_HPP_
