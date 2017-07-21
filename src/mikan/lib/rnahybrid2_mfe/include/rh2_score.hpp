#ifndef RH2_SCORE_HPP_
#define RH2_SCORE_HPP_

#include <vector>
#include <sstream>
#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"     // MKSiteScores
#include "rh2_option.hpp"        // RH2Options
#include "rh2_seed_site.hpp"     // RH2SeedSites
#include "hybrid_core.hpp"       // RH2WorkSpace, RH2RetValues

namespace rh2mfe {

//
// Store MFE scores
//
class RH2SiteScores : public mikan::MKSiteScores {
public:
    // Define methods
    RH2SiteScores(int pMaxMRNALen, int pMaxMiRNALen, std::string &pSeedDef) :
            mMaxMRNALen(pMaxMRNALen), mMaxMiRNALen(pMaxMiRNALen), mRHCore(pMaxMRNALen, pMaxMiRNALen, pSeedDef) {}

    void set_score(int i, float val) { mMFEScores[i] = val; };

    float const &get_score(int i) const { return mMFEScores[i]; }

    float const &get_norm_score(int i) const { return mNormScores[i]; }

    int get_hit_start(int pIdx) { return mRHRetVals[pIdx].mTargetPos0; }

    int get_hit_length(int pIdx) { return mRHRetVals[pIdx].mHitLen; }

    // Method prototype
    void init_from_args(RH2Options &opts) ;

    void clear_scores();

    int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites);

    void calc_normalized_score(int pIdx, int pTargetLen, int pQueryLen);

    void write_alignment(int pIdx, bool align_only = true);

private:
    int mMaxMRNALen;
    int mMaxMiRNALen;
    seqan::String<float> mMFEScores;
    seqan::String<float> mNormScores;
    rh2::RH2WorkSpace mRHCore;
    std::vector<rh2::RH2RetValues> mRHRetVals;
    seqan::CharString mOverlapDef;

private:
    void create_rh_seq(mikan::TRNAStr const &pRNASeq, std::vector<char> &pRHSeq);

    void write_seq_info(mikan::TRNAStr &pMiSeq, mikan::TRNAStr &pMRNASeq, std::vector<char> &pRhMiRNASeq,
                        std::vector<char> &pRhMRNASeq);
};

//
// Total MFE scores
//
class RH2TotalScores {
public:
    // Define methods
    RH2TotalScores() {}

    const seqan::String<float> &get_scores() { return mTotalScores; }

    const seqan::String<float> &get_norm_scores() { return mTotalNormScores; }

    const seqan::String<int> &get_mrna_pos() { return mMRNAPos; }

    const seqan::String<int> &get_site_num() { return mSiteNum; }

    // Method prototypes
    void clear_scores();

    int calc_scores(RH2SeedSites &pSeedSites, RH2SiteScores &pMFEScores,
                    const seqan::String<unsigned> &pSortedSites);

private:
    seqan::String<float> mTotalScores;
    seqan::String<float> mTotalNormScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace rh2mfe

#endif /* RH2_SCORE_HPP_ */
