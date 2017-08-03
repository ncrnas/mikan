#include <math.h>                // roundf
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_rna_score.hpp"     // TS5TotalScores

using namespace seqan;

namespace ts5cs {

//
// TS5TotalScores methods
//
void TS5TotalScores::clear_scores() {
    clear(mTotalScores);
    clear(mMRNAPos);
    clear(mSiteNum);
}

int TS5TotalScores::calc_scores(
        TS5SeedSites &pSeedSites,
        mikan::TRNASet const &,
        mikan::MKRMAWithSites &pRNAWithSites,
        TS5SiteScores &pContextScores) {

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    resize(mTotalScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs));
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs));

    float score, totalScore;
    unsigned siteCount;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        score = totalScore = 0;
        siteCount = 0;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            score = pContextScores.get_score(rnaSitePosMap[i][j]);
            totalScore += score;
            ++siteCount;
        }

        if (siteCount == 0) {
            continue;
        }

        mTotalScores[i] = totalScore;
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = siteCount;
    }

    return 0;

}

} // namespace ts5cs
