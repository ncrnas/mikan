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
    clear(mLogMaxEnScores);
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
    resize(mLogMaxEnScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mEffectiveRNAs, length(pRNAWithSites.mEffectiveRNAs), false);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);


    float score, maxScore, totalScore;
    float scoreEn, maxScoreEn, totalScoreEn;
    unsigned siteCount, maxIdx, maxIdxEn;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        score = scoreEn = totalScore = totalScoreEn = 0;
        maxScore = maxScoreEn = -FLT_MAX;
        siteCount = 0;
        maxIdx = maxScoreEn = 0;

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            score = pSiteScores.get_align_score(rnaSitePosMap[i][j]);
            totalScore += score;
            if (score > maxScore) {
                maxScore = score;
                maxIdx = j;
            }

            scoreEn = pSiteScores.get_energy_score(rnaSitePosMap[i][j]);
            totalScoreEn += scoreEn;
            if ((-1.0 * scoreEn) > maxScoreEn) {
                maxScoreEn = (-1.0 * scoreEn);
                maxIdxEn = j;
            }

            ++siteCount;
        }

        if (siteCount == 0) {
            continue;
        }

        mTotalAlignScores[i] = totalScore;
        mTotalEnScores[i] = totalScoreEn;
        mLogMaxAlignScores[i] = maxScore + std::log(totalScore);
        mLogMaxEnScores[i] = -1.0 * (maxScoreEn + std::log(-1.0 * totalScoreEn));
        mEffectiveRNAs[i] = true;
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
    }

    return 0;

}

} // namespace mr3as
