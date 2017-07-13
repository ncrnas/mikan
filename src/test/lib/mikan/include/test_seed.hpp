#ifndef MIKAN_TEST_SEED_HPP_
#define MIKAN_TEST_SEED_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"

template<class TSeedSeqs, class TTestIO>
class TestSeed : public TTestIO {
protected:

    void test_seed(const char *rnastr, int idx, const char *seq_type, bool effective, unsigned mmpos) {
        seqan::RnaString seedseq = rnastr;

        reverseComplement(seedseq);
        comp_two_rnas(mSeedSeqs.get_seed_seq(idx), seedseq);
        EXPECT_STREQ(seq_type, seqan::toCString((seqan::CharString) mSeedSeqs.get_seed_type(idx)));
        EXPECT_EQ(effective, mSeedSeqs.mEffectiveSeeds[idx]);
        EXPECT_EQ(mmpos, mSeedSeqs.get_mismatched_pos(idx));

//        seqan::RnaString seedseq2 = mSeedSeqs.get_seed_seq(idx);
//        reverseComplement(seedseq2);
//        std::cout << "SS Z" << seqan::toCString((seqan::CharString)seedseq2) << "Z, ";
//        std::cout << idx <<  ", ";
//        std::cout << "Z" << seqan::toCString((seqan::CharString)mSeedSeqs.get_seed_type(idx)) <<  "Z, ";
//        std::cout << mSeedSeqs.mEffectiveSeeds[idx] <<  ", ";
//        std::cout << mSeedSeqs.get_mismatched_pos(idx) <<  ");";
//        std::cout << std::endl;
    }

    void test_seed2(const char *rnastr, int idx, const char *seq_type, bool effective) {
        seqan::RnaString seedseq = rnastr;

        reverseComplement(seedseq);
        comp_two_rnas(mSeedSeqs.get_seed_seq(idx), seedseq);
        EXPECT_STREQ(seq_type, seqan::toCString((seqan::CharString) mSeedSeqs.get_seed_type(idx)));
        EXPECT_EQ(effective, mSeedSeqs.mEffectiveSeeds[idx]);
    }

    void test_seed3(const char *rnastr) {
        seqan::RnaString seedseq = rnastr;

        reverseComplement(seedseq);
        comp_two_rnas(mSeedSeqs.get_seed_seq(0), seedseq);
    }

    TSeedSeqs mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedDef;
    seqan::CharString mSeedDef1;
    seqan::CharString mOverlapDef;
};

#endif //MIKAN_TEST_SEED_HPP_
