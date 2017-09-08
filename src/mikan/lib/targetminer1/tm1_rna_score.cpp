#include <math.h>                // roundf
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_rna_score.hpp"     // TM1RNAScores

using namespace seqan;

namespace tm1p {

//
// TM1TotalScores methods
//
void TM1RNAScores::clear_scores() {
    mikan::MKRNAScores::clear_scores();

    mMRNAFeatures.clear_features();
    mMRNAInput.clear_scores();

    clear(mPredictions);
}

int TM1RNAScores::calc_scores(
        mikan::MKSeedSites &pSeedSites,
        mikan::TRNASet const &,
        mikan::MKRMAWithSites &pRNAWithSites,
        TM1SiteScores &mSiteScores) {

    int retVal;

    retVal = mMRNAFeatures.add_features(pSeedSites, pRNAWithSites, mSiteScores);
    if (retVal != 0) {
        return 1;
    }

    retVal = mMRNAInput.classify(mMRNAFeatures.get_scaled_feature());
    if (retVal != 0) {
        return 1;
    }

    const seqan::String<unsigned> &siteCounts = mMRNAFeatures.get_site_counts();
    const seqan::String<float> &scores = mMRNAInput.get_scores();
    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();

    resize(mRNAScores, length(siteCounts), 0.0);
    resize(mPredictions, length(siteCounts), 0);
    resize(mSiteNum, length(siteCounts), 0);
    resize(mMRNAPos, length(siteCounts), 0);
    resize(mEffectiveRNAs, length(siteCounts), false);

    for (unsigned i = 0; i < length(siteCounts); ++i) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        mRNAScores[i] = roundf(scores[i] * 1000.0) / 1000.0;
        if (mRNAScores[i] > 0) {
            mPredictions[i] = 1;
        } else {
            mPredictions[i] = -1;
        }
        mSiteNum[i] = siteCounts[i];
        mMRNAPos[i] = uniqRNAPosSet[i];
        mEffectiveRNAs[i] = true;
    }

    return 0;
}

} // namespace tm1p
