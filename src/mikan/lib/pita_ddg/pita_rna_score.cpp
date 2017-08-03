#include "mk_typedef.hpp"       // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "pita_site_score.hpp"  // PITASiteScores
#include "pita_rna_score.hpp"   // PITARNAScores

using namespace seqan;

namespace ptddg {

//
// PITARNAScores methods
//

int PITARNAScores::calc_scores(
        mikan::MKSeedSites &pSeedSites,
        mikan::TRNASet const &,
        mikan::MKRMAWithSites &pRNAWithSites,
        PITASiteScores &pSiteScores) {

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    resize(mRNAScores, length(pRNAWithSites.mEffectiveRNAs), 0.0);
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs));

    float score, max_score, exp_diff, total_score;
    unsigned site_count, max_idx;
    mikan::TSitePosSet sitePos;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        score = 0;
        max_score = -FLT_MAX;
        site_count = 0;
        max_idx = 0;

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }
            appendValue(sitePos, rnaSitePosMap[i][j]);
            score = -1.0 * pSiteScores.get_score(rnaSitePosMap[i][j]);
            if (score > max_score) {
                max_score = score;
                max_idx = j;
            }
            ++site_count;
        }

        if (length(sitePos) == 0) {
            continue;
        }

        total_score = 1.0;
        for (unsigned j = 0; j < length(sitePos); ++j) {
            if (j == max_idx) {
                continue;
            }
            score = -1.0 * pSiteScores.get_score(sitePos[j]);
            exp_diff = score - max_score;
            if (exp_diff > MIN_EXP_DIFF) {
                total_score += std::exp(exp_diff);
            }
        }

        mRNAScores[i] = -1.0 * (max_score + std::log(total_score));
        mMRNAPos[i] = uniqRNAPosSet[i];
        mSiteNum[i] = site_count;

        clear(sitePos);
    }

    return 0;

}

} // namespace ptddg
