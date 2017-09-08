#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_site_score.hpp"       // TM1SiteScores

using namespace seqan;

namespace tm1p {

//
// TSSVMSiteScores methods
//
void TM1SiteScores::clear_scores() {
    mikan::MKSiteScores::clear_scores();

    mSiteFeatures.clear_features();
}

int TM1SiteScores::calc_scores(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites &mRNAWithSites) {

    int retVal;

    // Generate site features
    retVal = mSiteFeatures.add_features(pMiRNASeq, pMRNASeqs, pSeedSites, mRNAWithSites);
    if (retVal != 0) {
        return 1;
    }

    mEffectiveSites = pSeedSites.mEffectiveSites;

    return 0;
}

} // namespace tm1p
