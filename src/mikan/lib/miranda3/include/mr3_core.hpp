#ifndef MR3_CORE_HPP_
#define MR3_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mr3_option.hpp"         // MR3Options
#include "mr3_score.hpp"          // MR3GGDScores, MR3TotalScores
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "mr3_site_cluster.hpp"   // MR3Overlap, MR3SortedSitePos

namespace mr3as {

int MR3CoreMain(int argc, char const **argv);

//
// MR3 score process core
//
class MR3Core {
public:
    // Constant values
    static const unsigned INDEXED_SEQ_LEN = mikan::SEEDLEN;
    static const unsigned OVERLAP_LEN = 6;

    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalSiteScore;
    bool mExecFilterOverlap;
    bool mExecSortSites;
    bool mExecSumScores;
    bool mOutputSiteScore;
    bool mOutputTotalScore;
    bool mOutputAlign;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileTotal;
    int mMinSeedLen;
    int mMaxSeedLen;
    float mMinAlignScore;
    float mMaxEnergy;

    mikan::TCharSet mSeedTypeDef;

public:
    // Define methods
    MR3Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecCalSiteScore(true), mExecFilterOverlap(true),
            mExecSortSites(true), mExecSumScores(true), mOutputSiteScore(true), mOutputTotalScore(true),
            mOutputAlign(true), mMinSeedLen(6), mMaxSeedLen(8), mMinAlignScore(120.0), mMaxEnergy(1.0),
            mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds),
            mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
            mSiteScores(pOpts) {
        init_from_args(pOpts);
    }

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

    MR3SeedSites mSeedSites;
    MR3SiteScores mSiteScores;
    mikan::MKRMAWithSites mRNAWithSites;
    MR3SiteFilter mSiteFilter;
    MR3TotalScores mTotalScores;

private:
    int write_site_score(seqan::CharString const &pMiRNAId);

    int write_total_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace mr3as

#endif /* MR3_CORE_HPP_ */
