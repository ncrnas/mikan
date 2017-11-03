#ifndef TSSVM_CORE_HPP_
#define TSSVM_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"          // MKSequences
#include "mk_option.hpp"            // MKOptions
#include "mk_core_tmpl.hpp"         // MKCoreTmpl
#include "tssvm_option.hpp"         // TSSVMOptions
#include "tssvm_seed_site.hpp"      // TSSVMSeedSites
#include "tssvm_site_filter.hpp"    // TSSVMSiteFilter
#include "tssvm_site_score.hpp"     // TSSVMSiteScores
#include "tssvm_rna_score.hpp"      // TSSVMRNAScores

namespace tssvm {

//
//Two-step SVM score process core
//
typedef mikan::MKCoreTmpl<TSSVMSeedSeqs, TSSVMSeedSites, TSSVMSiteScores, TSSVMSiteFilter, TSSVMRNAScores> TSSVMCoreBase;

class TSSVMCore : public TSSVMCoreBase {
public:
    // Define methods
    TSSVMCore(mikan::MKOptions const &pOpts,
              mikan::TCharSet const &pMiRNAIds,
              mikan::TRNASet const &pMiRNASeqs,
              mikan::TCharSet const &pMRNAIds,
              mikan::TRNASet const &pMRNASeqs,
              mikan::TIndexQGram &pRNAIdx,
              mikan::TFinder &pFinder) :
            TSSVMCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mFindSeedSites = true;
        mFilterSites = true;
        mCalcSiteScore = true;
        mClusterSites1 = true;
        mFilterSiteScores = false;
        mClusterSites2 = false;
        mSelectTopSites = false;
        mClusterSites3 = false;
        mCalcRNAScore = true;

    }

private:
    virtual void prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx);

    virtual void prepare_rna_output(mikan::TCharStr const &pMiRNAId);

    virtual int write_alignment(mikan::TCharStr const &pMiRNAId);

};

} // namespace tssvm

#endif /* TSSVM_CORE_HPP_ */
