#ifndef TSSVM_RNA_SCORE_HPP_
#define TSSVM_RNA_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_sites.hpp"         // MKRMAWithSites
#include "mk_rna_score.hpp"         // MKRNAScores
#include "tssvm_site_score.hpp"     // TSSVMSiteScores
#include "tssvm_mrna_feature.hpp"   // TSSVMRNARawFeatures
#include "tssvm_mrna_svm.hpp"       // TSSVMRNAInputVector

namespace tssvm {

//
// RNA level scores
//
class TSSVMRNAScores : public mikan::MKRNAScores {
public:
    // Define methods
    TSSVMRNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts) {}

    // Method prototypes
    void clear_scores();

    int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKRMAWithSites &pRNAWithSites, TSSVMSiteScores &pSitecores);

private:
    TSSVMRNARawFeatures mRnaFeatures;
    TSSVMRNAInputVector mRnaInput;

};

} // namespace tssvm

#endif /* TSSVM_RNA_SCORE_HPP_ */
