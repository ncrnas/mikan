#ifndef TS5_CORE_HPP_
#define TS5_CORE_HPP_

#include "mk_typedef.hpp"        // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "ts5_feature.hpp"       // TS5RawFeatures
#include "ts5_option.hpp"        // TS5CSOptions
#include "ts5_score.hpp"         // TS5ContextScores, TS5TotalScores
#include "ts5_seed_site.hpp"     // TS5SeedSites

namespace ts5cs {

int TS5CoreMain(int argc, char const **argv);

//
// TargetScan context score process core
//
class TS5Core {
public:
    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecGetRawFeat;
    bool mExecCalcContexScore;
    bool mExecSumScores;
    bool mOutputContexScore;
    bool mOutputTotalScore;
    bool mOutputAlign;
    seqan::CharString mOFileContext;
    seqan::CharString mOFileTotal;

public:
    // Define methods
    TS5Core(mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs, 
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecGetRawFeat(true), mExecCalcContexScore(true),
            mExecSumScores(true), mOutputContexScore(true), mOutputTotalScore(true),
            mOutputAlign(true), mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds),
            mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder) {}

    // Method prototypes
    void init_from_args(TS5CSOptions &opts);

    int open_output_file();

    int calculate_all_scores();

    int calculate_mirna_scores(unsigned pIdx);

private:
    mikan::TCharSet const &mMiRNAIds;
    mikan::TRNASet const &mMiRNASeqs;
    mikan::TCharSet const &mMRNAIds;
    mikan::TRNASet const &mMRNASeqs;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

    TS5SeedSites mSeedSites;
    TS5RawFeatures mRawFeatures;
    TS5ContextScores mCsScores;
    TS5TotalScores mTotalScore;

private:
    int write_context_score(seqan::CharString const &pMiRNAId);

    int write_total_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace ts5cs

#endif /* TS5_CORE_HPP_ */
