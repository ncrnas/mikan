#ifndef MK_CORE_BASE_HPP_
#define MK_CORE_BASE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_option.hpp"          // MKOptions
#include "mk_seed_site.hpp"       // MKSeedSites
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_score.hpp"       // MKRNAScores

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
    bool mPrintSiteHeader;
    bool mPrintRNAheader;

    mikan::TCharStr mOFileSite;
    mikan::TCharStr mOFileRNA;

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
            mPrintSiteHeader(true),
            mPrintRNAheader(true),
            mOpts(pOpts),
            mMiRNAIds(pMiRNAIds), mMiRNASeqs(pMiRNASeqs),
            mMRNAIds(pMRNAIds), mMRNASeqs(pMRNASeqs),
            mRNAIdx(pRNAIdx), mFinder(pFinder) {

        init_from_args(pOpts);

    }

    // Method prototypes
    virtual void init_from_args(mikan::MKOptions const &opts);

    int open_output_file();

    void close_output_file();

    int calculate_all_scores();

    int calculate_mirna_scores(unsigned pIdx);

    virtual int find_seed_sites(unsigned pIdx) = 0;

    virtual int calc_site_scores(unsigned pIdx) = 0;

    virtual int calc_rna_scores(unsigned pIdx) = 0;

    virtual int output_results(unsigned pIdx) = 0;

    virtual void clear_all() = 0;

    virtual mikan::MKSeedSites &get_seed_sites() = 0;

    virtual mikan::MKSiteScores &get_site_scores() = 0;

    virtual mikan::MKRNAScores &get_rna_scores() = 0;

protected:
    mikan::MKOptions const &mOpts;
    mikan::TCharSet const &mMiRNAIds;
    mikan::TRNASet const &mMiRNASeqs;
    mikan::TCharSet const &mMRNAIds;
    mikan::TRNASet const &mMRNASeqs;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;

    std::ofstream mOFile1;
    std::ofstream mOFile2;

    int write_site_score(mikan::TCharStr const &pMiRNAId, mikan::MKSeedSites &pSeedSites,
                         mikan::MKRMAWithSites &pRNAWithSites);

    virtual void prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) = 0;

    void write_site_score_tab(std::string &pMiRNAName, std::string &pMRNAName, unsigned pStartPos, unsigned pEndPos,
                              std::string &pSeedType, std::string &pScore1Name, std::string &pScore1,
                              std::string &pScore2Name, std::string &pScore2);

    void write_site_score_gff(std::string &pHeader, std::string &pSrc, std::string &pMiRNAName, std::string &pMRNAName,
                              unsigned pStartPos, unsigned pEndPos, std::string &pSeedType,
                              std::string &pScore1, std::string &pScore2);

    int write_rna_score(mikan::TCharStr const &pMiRNAId);

    virtual void prepare_rna_output(mikan::TCharStr const &pMiRNAId) = 0;

    void write_rna_score_tab(std::string &pMiRNAName, std::string &pMRNAName, unsigned pSiteNum,
                             std::string &pScore1Name, std::string &pScore1,
                             std::string &pScore2Name, std::string &pScore2);

    void write_rna_score_gff(std::string &pHeader, std::string &pSrc, std::string &pMiRNAName, std::string &pMRNAName,
                             unsigned pMRNALen, unsigned pSiteNum, std::string &pScore1, std::string &pScore2);

    virtual int write_alignment(mikan::TCharStr const &pMiRNAId) = 0;


};

} // namespace mikan

#endif /* MK_CORE_BASE_HPP_ */
