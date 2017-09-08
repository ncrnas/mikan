
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_rna_score.hpp"   // TSSVMRNAScores

using namespace seqan;

namespace tssvm {

//
// TSSVMRNAScores methods
//
void TSSVMRNAScores::clear_scores() {
    mikan::MKRNAScores::clear_scores();

    mRnaFeatures.clear_features();
    mRnaInput.clear_scores();
}

int TSSVMRNAScores::calc_scores(
        mikan::MKSeedSites &pSeedSites,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKRMAWithSites &pRNAWithSites,
        TSSVMSiteScores &pSitecores) {

    int retVal;
    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();

    resize(mEffectiveRNAs, length(pRNAWithSites.mEffectiveRNAs), true);
    resize(mRNAScores, length(pRNAWithSites.mEffectiveRNAs));
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs));
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs));

    retVal = mRnaFeatures.add_features(pSeedSites, pMRNASeqs, pRNAWithSites, pSitecores);
    if (retVal != 0) {
        return 1;
    }

    retVal = mRnaInput.classify(mRnaFeatures);
    if (retVal != 0) {
        return 1;
    }

    mRNAScores = mRnaInput.get_scores();
    mMRNAPos = uniqRNAPosSet;
    mSiteNum = mRnaFeatures.get_site_count();

    return 0;

}

} // namespace tssvm
