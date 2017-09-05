#include "mk_typedef.hpp"       // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mr3_rna_score.hpp"    // MR3RNAScores

using namespace seqan;

namespace mr3as {

//
// MR3TotalScores methods
//
void MR3RNAScores::clear_scores() {
    mikan::MKRNAScores::clear_scores();

    clear(mTotalEnScores);
    clear(mTotalAlignScores);
    clear(mLogMaxAlignScores);
    clear(mLogMinEnScores);
    clear(mMaxAlignScores);
    clear(mMinEnScores);
}

int MR3RNAScores::calc_scores(
        mikan::MKSeedSites &pSeedSites,
        mikan::TRNASet const &,
        mikan::MKRMAWithSites &pRNAWithSites,
        MR3SiteScores &pSiteScores) {

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    resize(mTotalAlignScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mTotalEnScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mLogMaxAlignScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mLogMinEnScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMaxAlignScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMinEnScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mEffectiveRNAs, length(pRNAWithSites.mEffectiveRNAs), false);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);

    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        unsigned siteCount = 0;
        float totalScore = 0;
        float totalScoreEn = 0;
        float maxScore = -FLT_MAX;
        float minScoreEn = FLT_MAX;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            float score = pSiteScores.get_align_score(rnaSitePosMap[i][j]);
            totalScore += score;
            if (score > maxScore) {
                maxScore = score;
            }

            float scoreEn = pSiteScores.get_energy_score(rnaSitePosMap[i][j]);
            totalScoreEn += scoreEn;
            if (scoreEn < minScoreEn) {
                minScoreEn = scoreEn;
            }

            ++siteCount;
        }

        if (siteCount == 0) {
            continue;
        }

        mTotalAlignScores[i] = totalScore;
        mTotalEnScores[i] = totalScoreEn;
        if (totalScore - maxScore > 0) {
            mLogMaxAlignScores[i] = maxScore + std::log(totalScore - maxScore);
        } else {
            mLogMaxAlignScores[i] = maxScore;
        }
        if (totalScoreEn - minScoreEn < 0) {
            mLogMinEnScores[i] = minScoreEn + (-1.0 * std::log(-1.0 * (totalScoreEn - minScoreEn)));
        } else {
            mLogMinEnScores[i] = minScoreEn;
        }
        mMaxAlignScores[i] = maxScore;
        mMinEnScores[i] = minScoreEn;
        mEffectiveRNAs[i] = true;
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
    }

    return 0;

}

} // namespace mr3as
