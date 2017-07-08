#include <math.h>                // roundf
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_score.hpp"         // TM1ClassifiedScores

using namespace seqan;

namespace tm1p {

//
// TM1TotalScores methods
//
void TM1ClassifiedScores::clear_scores() {
    clear(mScores);
    clear(mPredictions);
    clear(mSiteNum);
}

int TM1ClassifiedScores::calc_scores(
        const seqan::String<unsigned> &pSiteCoutns,
        const seqan::String<float> &pScores) {
    resize(mScores, length(pSiteCoutns));
    resize(mPredictions, length(pSiteCoutns));
    resize(mSiteNum, length(pSiteCoutns));

    for (unsigned i = 0; i < length(pSiteCoutns); ++i) {

        mScores[i] = roundf(pScores[i] * 1000.0) / 1000.0;
        if (pScores[i] > 0) {
            mPredictions[i] = 1;
        } else {
            mPredictions[i] = -1;
        }
        mSiteNum[i] = pSiteCoutns[i];
    }

    return 0;
}

} // namespace tm1p
