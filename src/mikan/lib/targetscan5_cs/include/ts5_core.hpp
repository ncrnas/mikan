#ifndef TS5_CORE_HPP_
#define TS5_CORE_HPP_

#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "ts5_feature.hpp"        // TS5RawFeatures
#include "ts5_option.hpp"         // TS5CSOptions
#include "ts5_site_score.hpp"     // TS5SiteScores
#include "ts5_seed_site.hpp"      // TS5SeedSites
#include "ts5_rna_score.hpp"      // TS5RNAScores

namespace ts5cs {

int TS5CoreMain(int argc, char const **argv);

//
// TargetScan context score process core
//
class TS5Core {
public:
    // Declare variables
    bool mExecSearchSeedSites;
    bool mExecCalcSiteScore;
    bool mExecSumScores;
    bool mOutputContexScore;
    bool mOutputTotalScore;
    bool mOutputAlign;

    seqan::CharString mOFileSite;
    seqan::CharString mOFileRNA;

public:
    // Define methods
    TS5Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            mExecSearchSeedSites(true), mExecCalcSiteScore(true),
            mExecSumScores(true), mOutputContexScore(true), mOutputTotalScore(true),
            mOutputAlign(true), mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds),
            mMRNASeqs(pMRNASeqs), mSeedSeqs(pOpts), mSeedSites(pRNAIdx, pFinder, pMRNASeqs), mSiteScores(pOpts),
            mRNAWithSites(pOpts), mRNAScores(pOpts) {

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

    TS5SeedSeqs mSeedSeqs;
    TS5SeedSites mSeedSites;
    TS5SiteScores mSiteScores;
    mikan::MKRMAWithSites mRNAWithSites;
    TS5RNAScores mRNAScores;

private:
    int write_site_score(seqan::CharString const &pMiRNAId);

    int write_rna_score(seqan::CharString const &pMiRNAId);

    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace ts5cs

#endif /* TS5_CORE_HPP_ */
