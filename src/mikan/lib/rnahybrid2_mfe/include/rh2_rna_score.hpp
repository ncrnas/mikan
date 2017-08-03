#ifndef RH2_RNA_SCORE_HPP_
#define RH2_RNA_SCORE_HPP_

#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "rh2_site_score.hpp"     // RH2SiteScores

namespace rh2mfe {

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
                    mikan::MKRMAWithSites &pRNAWithSites);

private:
    seqan::String<float> mTotalScores;
    seqan::String<float> mTotalNormScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace rh2mfe

#endif /* RH2_RNA_SCORE_HPP_ */
