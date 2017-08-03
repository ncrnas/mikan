#ifndef TM1_CORE_HPP_
#define TM1_CORE_HPP_

#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "mk_site_score.hpp"     // MKSiteScores
#include "mk_option.hpp"         // MKOptions
#include "tm1_site_filter.hpp"   // TM1SiteFilter
#include "tm1_mrna_feature.hpp"  // TM1MRNAFeatures
#include "tm1_mrna_svm.hpp"      // TM1MRNAModel, TM1MRNAInputVector
#include "tm1_option.hpp"        // TM1CSOptions
#include "tm1_rna_score.hpp"     // TM1ClassifiedScores
#include "tm1_seed_site.hpp"     // TM1SeedSites
#include "tm1_site_score.hpp"    // TM1SiteScores

namespace tm1p {

int TM1CoreMain(int argc, char const **argv);

//
// TargetScan context score process core
//
class TM1Core {
public:
    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalSiteScore;
    bool mExecCalcSiteScore;
    bool mExecFilterOverlap;
    bool mExecSortSites;
    bool mExecGetMRNAFeat;
    bool mExecRNAScore;
    bool mExecSumScores;
    bool mOutputSitePos;
    bool mOutputScore;
    bool mOutputAlign;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileScore;

    mikan::TCharSet mSeedTypeDef;

public:
    // Define methods
    TM1Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecCalSiteScore(false), mExecCalcSiteScore(true), mExecFilterOverlap(true),
            mExecSortSites(true),
            mExecGetMRNAFeat(true), mExecRNAScore(true), mExecSumScores(true), mOutputSitePos(true),
            mOutputScore(true), mOutputAlign(true), mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds),
            mMRNASeqs(pMRNASeqs), mSeedSites(pRNAIdx, pFinder, pMRNASeqs), mSiteScores(pOpts), mSiteFilter(pOpts),
            mRNAScores(pOpts) {
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

    TM1SeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    TM1SiteScores mSiteScores;
    TM1SiteFilter mSiteFilter;
    TM1ClassifiedScores mRNAScores;

private:
    int write_site_positions(seqan::CharString const &pMiRNAId);

    int write_scores(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tm1p

#endif /* TM1_CORE_HPP_ */
