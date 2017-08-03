#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_option.hpp"         // TSSVMOptions
#include "tssvm_site_score.hpp"     // TSSVMSiteScores

using namespace seqan;

namespace tssvm {

//
// TSSVMSiteScores methods
//
void TSSVMSiteScores::clear_scores() {
    mikan::MKSiteScores::clear_scores();

    mAlignSeqs.clear_alignments();
    mSiteFeatures.clear_features();
    mSiteInput.clear_scores();
}

int TSSVMSiteScores::calc_scores(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites const &) {

    int retVal;

    // Align sequences
    retVal = mAlignSeqs.align_seq(pMiRNASeq, pMRNASeqs, pSeedSites);
    if (retVal != 0) {
        return 1;
    }

    // Generate site features
    retVal = mSiteFeatures.add_features(pMiRNASeq, pMRNASeqs, pSeedSites, mAlignSeqs);
    if (retVal != 0) {
        return 1;
    }

    // Calculate site SVM scores
    retVal = mSiteInput.classify(mSiteFeatures);
    if (retVal != 0) {
        return 1;
    }

    mEffectiveSites = pSeedSites.mEffectiveSites;

    return 0;
}

} // namespace tssvm
