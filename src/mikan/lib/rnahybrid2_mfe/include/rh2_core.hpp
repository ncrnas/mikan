#ifndef RH2_CORE_HPP_
#define RH2_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "rh2_option.hpp"        // RH2Options
#include "rh2_score.hpp"         // RH2MFEScores, RH2TotalScores
#include "rh2_seed_site.hpp"     // RH2SeedSites
#include "rh2_site_cluster.hpp"  // RH2Overlap, RH2TopNScore, RH2SortedSitePos

namespace rh2mfe {

int RH2CoreMain(int argc, char const **argv);

//
// RNAhybrid MFE score process core
//
class RH2Core {
public:
    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalMFEScore;
    bool mExecFilterOverlap;
    bool mExecFilterSiteNum;
    bool mExecSortSites;
    bool mExecSumScores;
    bool mOutputMFEScore;
    bool mOutputTotalScore;
    bool mOutputAlign;
    seqan::CharString mOFileMFE;
    seqan::CharString mOFileTotal;

    seqan::CharString mSeedDef;
    seqan::CharString mOverlapDef;
    int mMaxHits;

public:
    // Define methods
    RH2Core(mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, int pMRNAMaxLen, int pMiRNAMaxLen,
            std::string &pSeedDef) :
            mExecSearchSeedSites(true), mExecCalMFEScore(true), mExecFilterOverlap(true),
            mExecFilterSiteNum(true), mExecSortSites(true), mExecSumScores(true), mOutputMFEScore(true),
            mOutputTotalScore(true), mOutputAlign(true), mMaxHits(0), mMiRNAIds(pMiRNAIds),
            mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs),
            mSeedSites(pRNAIdx, pFinder, pMRNASeqs), mMfeScores(pMRNAMaxLen, pMiRNAMaxLen, pSeedDef) {}

    // Method prototypes
    void init_from_args(RH2Options &opts);

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

    RH2SeedSites mSeedSites;
    RH2MFEScores mMfeScores;
    RH2Overlap mOverlappedSites;
    RH2TopNScore mTopScoredSites;
    RH2SortedSitePos mSortedSites;
    RH2TotalScores mTotalScores;

private:
    int write_mfe_score(seqan::CharString const &pMiRNAId);

    int write_total_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace rh2mfe

#endif /* RH2_CORE_HPP_ */
