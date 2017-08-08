#ifndef RH2_CORE_HPP_
#define RH2_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "mk_option.hpp"         // MKOptions
#include "rh2_option.hpp"        // RH2Options
#include "rh2_site_score.hpp"    // RH2SiteScores
#include "rh2_seed_site.hpp"     // RH2SeedSites
#include "rh2_site_filter.hpp"   // RH2SiteFilter, RH2TopNSites
#include "rh2_rna_score.hpp"     // RH2RNAScores

namespace rh2mfe {

int RH2CoreMain(int argc, char const **argv);

//
// RNAhybrid MFE score process core
//
class RH2Core {
public:
    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalSiteScore;
    bool mExecFilterOverlap;
    bool mExecFilterSiteNum;
    bool mExecSortSites;
    bool mExecSumScores;
    bool mOutputMFEScore;
    bool mOutputTotalScore;
    bool mOutputAlign;

    seqan::CharString mOFileSite;
    seqan::CharString mOFileRNA;

    seqan::CharString mOverlapDef;
    int mMaxHits;

public:
    // Define methods
    RH2Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecCalSiteScore(true), mExecFilterOverlap(true),
            mExecFilterSiteNum(true), mExecSortSites(true), mExecSumScores(true), mOutputMFEScore(true),
            mOutputTotalScore(true), mOutputAlign(true), mMaxHits(0), mMiRNAIds(pMiRNAIds),
            mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs),
            mSeedSeqs(pOpts), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),  mRNAWithSites(pOpts), mSiteScores(pOpts),
            mSiteFilter(pOpts), mTopNSites(pOpts), mRNAScores(pOpts) {

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

    RH2SeedSeqs mSeedSeqs;
    RH2SeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    RH2SiteScores mSiteScores;
    RH2SiteFilter mSiteFilter;
    mikan::MKTopNSites mTopNSites;
    RH2RNAScores mRNAScores;

private:
    int write_site_score(seqan::CharString const &pMiRNAId);

    int write_rna_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace rh2mfe

#endif /* RH2_CORE_HPP_ */
