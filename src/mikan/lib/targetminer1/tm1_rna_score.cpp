#include <math.h>                // roundf
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_rna_score.hpp"     // TM1ClassifiedScores

using namespace seqan;

namespace tm1p {

//
// TM1TotalScores methods
//
void TM1ClassifiedScores::clear_scores() {
    mikan::MKRNAScores::clear_scores();

    mMRNAFeatures.clear_features();
    mMRNAInput.clear_scores();

    clear(mPredictions);
}

int TM1ClassifiedScores::calc_scores(
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

    resize(mRNAScores, length(siteCounts));
    resize(mPredictions, length(siteCounts));
    resize(mSiteNum, length(siteCounts));

    for (unsigned i = 0; i < length(siteCounts); ++i) {

        mRNAScores[i] = roundf(scores[i] * 1000.0) / 1000.0;
        if (mRNAScores[i] > 0) {
            mPredictions[i] = 1;
        } else {
            mPredictions[i] = -1;
        }
        mSiteNum[i] = siteCounts[i];
    }

    return 0;
}

} // namespace tm1p
