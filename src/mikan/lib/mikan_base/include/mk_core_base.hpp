#ifndef MK_CORE_BASE_HPP_
#define MK_CORE_BASE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_option.hpp"          // MKOptions

namespace mikan {

//
// Base class of MK process core base
//
class MKCoreBase {
public:
    // Declare variables
    bool mFindSeedSites;
    bool mFilterSites;
    bool mCalcSiteScore;
    bool mClusterSites1;
    bool mFilterSiteScores;
    bool mClusterSites2;
    bool mSelectTopSites;
    bool mClusterSites3;
    bool mCalcRNAScore;
    bool mOutputSite;
    bool mOutputRNA;
    bool mOutputAlign;

    seqan::CharString mOFileSite;
    seqan::CharString mOFileRNA;

    // Define methods
    MKCoreBase(mikan::MKOptions const &pOpts,
               mikan::TCharSet const &pMiRNAIds,
               mikan::TRNASet const &pMiRNASeqs,
               mikan::TCharSet const &pMRNAIds,
               mikan::TRNASet const &pMRNASeqs,
               mikan::TIndexQGram &pRNAIdx,
               mikan::TFinder &pFinder) :
            mFindSeedSites(true),
            mFilterSites(true),
            mCalcSiteScore(true),
            mClusterSites1(true),
            mFilterSiteScores(true),
            mClusterSites2(true),
            mSelectTopSites(true),
            mClusterSites3(true),
            mCalcRNAScore(true),
            mOutputSite(true),
            mOutputRNA(true),
            mOutputAlign(true),
            mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs), mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs),
            mRNAIdx(pRNAIdx), mFinder(pFinder) {

        init_from_args(pOpts);

    }

    // Method prototypes
    void init_from_args(mikan::MKOptions const &opts);

    int open_output_file();

    int calculate_all_scores();

    int calculate_mirna_scores(unsigned pIdx);

    virtual int find_seed_sites(unsigned pIdx) = 0;

    virtual int calc_site_scores(unsigned pIdx) = 0;

    virtual int calc_rna_scores(unsigned pIdx) = 0;

    virtual int output_results(unsigned pIdx) = 0;

    virtual void clear_all() = 0;

protected:
    mikan::TCharSet const &mMiRNAIds;
    mikan::TRNASet const &mMiRNASeqs;
    mikan::TCharSet const &mMRNAIds;
    mikan::TRNASet const &mMRNASeqs;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

};

} // namespace mikan

#endif /* MK_CORE_BASE_HPP_ */
