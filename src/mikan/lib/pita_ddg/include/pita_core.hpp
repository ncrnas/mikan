#ifndef PITA_CORE_HPP_
#define PITA_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_core_base.hpp"       // MKCoreBase
#include "pita_option.hpp"        // PITAOptions
#include "pita_site_score.hpp"    // PITAGGDScores
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_filter.hpp"   // PITASiteFilter
#include "pita_rna_score.hpp"     // PITARNAScores

namespace ptddg {

int PITACoreMain(int argc, char const **argv);

//
// PITA score process core
//
class PITACore : public mikan::MKCoreBase {
public:
    // Define methods
    PITACore(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
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
    PITASeedSeqs mSeedSeqs;
    PITASeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    PITASiteScores mSiteScores;
    PITASiteFilter mSiteFilter;
    mikan::MKTopNSites mTopNSites;
    PITARNAScores mRNAScores;

    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_rna_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace ptddg

#endif /* PITA_CORE_HPP_ */
