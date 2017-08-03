#ifndef PITA_CORE_HPP_
#define PITA_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "pita_option.hpp"        // PITAOptions
#include "pita_site_score.hpp"    // PITAGGDScores
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_filter.hpp"   // PITASiteFilter
#include "pita_rna_score.hpp"     // PITARNAScores

namespace ptddg {

int PITACoreMain(int argc, char const **argv);

//
// PITA score process core
//
class PITACore {
public:
    // Constant values
    static const unsigned INDEXED_SEQ_LEN = mikan::SEEDLEN;
    static const unsigned OVERLAP_LEN = 0;

    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalSiteScore;
    bool mExecFilterOverlap;
    bool mExecSortSites;
    bool mExecSumScores;
    bool mOutputDDGScore;
    bool mOutputTotalScore;
    bool mOutputAlign;
    seqan::CharString mOFileDDG;
    seqan::CharString mOFileTotal;
    int mMinSeedLen;
    int mMaxSeedLen;

    mikan::TCharSet mSeedTypeDef;

public:
    // Define methods
    PITACore(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
             mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
             mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecCalSiteScore(true), mExecFilterOverlap(true),
            mExecSortSites(true), mExecSumScores(true), mOutputDDGScore(true), mOutputTotalScore(true),
            mOutputAlign(true), mMinSeedLen(6), mMaxSeedLen(8),
            mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds),
            mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
            mSiteScores(pOpts), mSiteFilter(pOpts), mRNAScores(pOpts) {
        init_from_args(pOpts);
    }

    void set_backtrack(bool pBT) { mSiteScores.set_backtrack(pBT); }

    // Method prototypes
    void init_from_args(mikan::MKOptions const &opts);

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

    PITASeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    PITASiteScores mSiteScores;
    PITASiteFilter mSiteFilter;
    PITARNAScores mRNAScores;

private:
    int write_ddg_score(seqan::CharString const &pMiRNAId);

    int write_total_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace ptddg

#endif /* PITA_CORE_HPP_ */
