#include "mk_typedef.hpp"      // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "rh2_rna_score.hpp"   // RH2TotalScores

using namespace seqan;

namespace rh2mfe {

//
// RH2TotalScores methods
//
void RH2TotalScores::clear_scores() {
    clear(mMRNAPos);
    clear(mSiteNum);
    clear(mTotalScores);
    clear(mTotalNormScores);
}

int RH2TotalScores::calc_scores(
        RH2SeedSites &pSeedSites,
        RH2SiteScores &pMFEScores,
        mikan::MKRMAWithSites &pRNAWithSites) {

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    resize(mTotalScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mTotalNormScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs));
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);


    float totalScore, totalNormScore;
    unsigned siteCount;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        totalScore = 0;
        totalNormScore = 0;
        siteCount = 0;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            totalScore += pMFEScores.get_score(rnaSitePosMap[i][j]);
            totalNormScore += pMFEScores.get_norm_score(rnaSitePosMap[i][j]);
            ++siteCount;
        }

        if (siteCount == 0) {
            continue;
        }

        mTotalScores[i] = totalScore;
        mTotalNormScores[i] = totalNormScore;
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
    }

    return 0;

}

} // namespace rh2mfe
