#ifndef TSSVM_CORE_HPP_
#define TSSVM_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"          // MKSequences
#include "mk_option.hpp"            // MKOptions
#include "mk_core.hpp"              // MKCoreBase
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
class TSSVMCore : public mikan::MKCoreBase  {
public:
    // Define methods
    TSSVMCore(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
    mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
    mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
    MKCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
    mSeedSeqs(pOpts), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
    mRNAWithSites(pOpts), mSiteScores(pOpts), mSiteFilter(pOpts), mTopNSites(pOpts), mRNAScores(pOpts) {

        mFilterSiteScores = false;
        mSelectTopSites = false;
        mSelectTopRNAs = false;

    }

    virtual int find_seed_sites(unsigned pIdx);
    virtual int calc_site_scores(unsigned pIdx);
    virtual int ensemble_site_scores(unsigned pIdx);
    virtual int calc_rna_scores(unsigned pIdx);
    virtual int ensemble_rna_scores(unsigned pIdx);
    virtual int output_results(unsigned pIdx);
    virtual void clear_all();

private:
    TSSVMSeedSeqs mSeedSeqs;
    TSSVMSeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    TSSVMSiteScores mSiteScores;
    TSSVMSiteFilter mSiteFilter;
    mikan::MKTopNSites mTopNSites;
    TSSVMRNAScores mRNAScores;
    
    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_rna_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tssvm

#endif /* TSSVM_CORE_HPP_ */
