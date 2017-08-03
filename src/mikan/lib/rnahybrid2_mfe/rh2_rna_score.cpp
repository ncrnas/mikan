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
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs));
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);


    float rnaScore, normScore;
    unsigned siteCount;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        rnaScore = 0;
        normScore = 0;
        siteCount = 0;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            rnaScore += pMFEScores.get_score(rnaSitePosMap[i][j]);
            normScore += pMFEScores.get_norm_score(rnaSitePosMap[i][j]);
            ++siteCount;
        }

        if (siteCount == 0) {
            continue;
        }

        mRNAScores[i] = rnaScore;
        mNormScores[i] = normScore;
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
    }

    return 0;

}

} // namespace rh2mfe
