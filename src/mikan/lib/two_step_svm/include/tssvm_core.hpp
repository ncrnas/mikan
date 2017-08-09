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

int TSSVMCoreMain(int argc, char const **argv);

//
//Two-step SVM score process core
//
typedef mikan::MKCoreTmpl<TSSVMSeedSeqs, TSSVMSeedSites, TSSVMSiteScores, TSSVMSiteFilter, TSSVMRNAScores > TSSVMCoreBase;

class TSSVMCore : public TSSVMCoreBase  {
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

        mClusterSites2 = false;
        mFilterSiteScores = false;
        mSelectTopSites = false;
        mClusterSites3 = false;

    }

private:
    virtual int write_site_score(seqan::CharString const &pMiRNAId);
    virtual int write_rna_score(seqan::CharString const &pMiRNAId);
    virtual int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tssvm

#endif /* TSSVM_CORE_HPP_ */
