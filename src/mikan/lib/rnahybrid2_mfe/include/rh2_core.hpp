#ifndef RH2_CORE_HPP_
#define RH2_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "mk_option.hpp"         // MKOptions
#include "mk_core.hpp"           // MKCoreBase
#include "rh2_option.hpp"        // RH2Options
#include "rh2_site_score.hpp"    // RH2SiteScores
#include "rh2_seed_site.hpp"     // RH2SeedSites
#include "rh2_site_filter.hpp"   // RH2SiteFilter, RH2TopNSites
#include "rh2_rna_score.hpp"     // RH2RNAScores

namespace rh2mfe {

int RH2CoreMain(int argc, char const **argv);

//
// RNAhybrid MFE score process core
//
class RH2Core : public mikan::MKCoreBase  {
public:
    // Define methods
    RH2Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
    mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
    mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
    MKCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
    mSeedSeqs(pOpts), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),  mRNAWithSites(pOpts), mSiteScores(pOpts),
    mSiteFilter(pOpts), mTopNSites(pOpts), mRNAScores(pOpts) {

        mFilterSites = false;
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
    RH2SeedSeqs mSeedSeqs;
    RH2SeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    RH2SiteScores mSiteScores;
    RH2SiteFilter mSiteFilter;
    mikan::MKTopNSites mTopNSites;
    RH2RNAScores mRNAScores;

    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_rna_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace rh2mfe

#endif /* RH2_CORE_HPP_ */
