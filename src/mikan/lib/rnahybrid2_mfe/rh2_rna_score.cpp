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
    clear(mLogMinRNAScores);
    clear(mLogMaxNormScores);
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
    resize(mLogMinRNAScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mLogMaxNormScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mEffectiveRNAs, length(pRNAWithSites.mEffectiveRNAs), false);

    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        unsigned siteCount = 0;
        float totRnaScore = 0;
        float totNormScore = 0;
        float minRnaScore= FLT_MAX;
        float maxNormScore = -FLT_MAX;

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            float score = pMFEScores.get_score(rnaSitePosMap[i][j]);
            if (score < minRnaScore) {
                minRnaScore = score;
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
        mMinScores[i] = minRnaScore;
        mMaxNormScores[i] = maxNormScore;
        if (totRnaScore - minRnaScore < 0) {
            mLogMinRNAScores[i] = minRnaScore + (-1.0 * std::log(-1.0 * (totRnaScore - minRnaScore)));

        } else {
            mLogMinRNAScores[i] = minRnaScore;
        }
        if (totNormScore - maxNormScore > 0) {
            mLogMaxNormScores[i] = maxNormScore + std::log(totNormScore - maxNormScore);
        } else {
            mLogMaxNormScores[i] = maxNormScore;
        }
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
        mEffectiveRNAs[i] = true;
    }

    return 0;

}

} // namespace rh2mfe
