#include "mk_typedef.hpp"      // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "rh2_rna_score.hpp"   // RH2RNAScores

using namespace seqan;

namespace rh2mfe {

//
// RH2RNAScores methods
//
void RH2RNAScores::clear_scores() {
    mikan::MKRNAScores::clear_scores();

    clear(mNormScores);
    clear(mMinScores);
    clear(mMaxNormScores);
}

int RH2RNAScores::calc_scores(
        mikan::MKSeedSites &pSeedSites,
        mikan::TRNASet const &,
        mikan::MKRMAWithSites &pRNAWithSites,
        RH2SiteScores &pMFEScores) {

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    resize(mRNAScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mNormScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMinScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMaxNormScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mEffectiveRNAs, length(pRNAWithSites.mEffectiveRNAs), false);

    float totRnaScore, totNormScore;
    unsigned siteCount;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        totRnaScore = 0;
        totNormScore = 0;
        siteCount = 0;
        float maxScore, maxNormScore;
        maxScore = maxNormScore = -FLT_MAX;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            float score = pMFEScores.get_score(rnaSitePosMap[i][j]);
            if ((-1.0 * score) > maxScore) {
                maxScore = -1.0 * score;
            }

            float normScore = pMFEScores.get_norm_score(rnaSitePosMap[i][j]);
            if (normScore > maxNormScore) {
                maxNormScore = normScore;
            }

            totRnaScore += score;
            totNormScore += normScore;
            ++siteCount;
        }

        if (siteCount == 0) {
            continue;
        }

        mRNAScores[i] = totRnaScore;
        mNormScores[i] = totNormScore;
        mMinScores[i] = -1.0 * maxScore;
        mMaxNormScores[i] = maxNormScore;
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
        mEffectiveRNAs[i] = true;
    }

    return 0;

}

} // namespace rh2mfe
