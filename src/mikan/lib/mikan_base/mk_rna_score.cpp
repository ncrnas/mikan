#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_score.hpp"      // MKRNAScores

using namespace seqan;

namespace mikan {

//
// MKSiteScores methods
//
void MKRNAScores::clear_scores() {
    clear(mEffectiveRNAs);
    clear(mRNAScores);
    clear(mMRNAPos);
    clear(mSiteNum);

}

int MKRNAScores::calc_scores(
        mikan::MKSeedSites &pSeedSites,
        mikan::TRNASet const &,
        mikan::MKRMAWithSites &pRNAWithSites,
        mikan::MKSiteScores &pSiteScores) {

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    resize(mEffectiveRNAs, length(pRNAWithSites.mEffectiveRNAs), false);
    resize(mRNAScores, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);


    float rnaScore;
    unsigned siteCount;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        rnaScore = 0;
        siteCount = 0;

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            rnaScore += pSiteScores.get_score(rnaSitePosMap[i][j]);
            ++siteCount;
        }

        if (siteCount == 0) {
            continue;
        }

        mRNAScores[i] = rnaScore;
        mEffectiveRNAs[i] = true;
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
    }

    return 0;
}

} // namespace mikan
