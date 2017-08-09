#ifndef MR3_CORE_HPP_
#define MR3_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_core_base.hpp"       // MKCoreBase
#include "mr3_option.hpp"         // MR3Options
#include "mr3_site_score.hpp"     // MR3GGDScores
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "mr3_site_filter.hpp"    // MR3SiteFilter
#include "mr3_rna_score.hpp"      // MR3RNAScores

namespace mr3as {

int MR3CoreMain(int argc, char const **argv);

//
// MR3 score process core
//
class MR3Core : public mikan::MKCoreBase {
public:
    // Define methods
    MR3Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
            MKCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mSeedSeqs(pOpts), mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
            mSiteScores(pOpts), mRNAWithSites(pOpts), mSiteFilter(pOpts), mTopNSites(pOpts), mRNAScores(pOpts) {


        mFilterSites = false;
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
    MR3SeedSeqs mSeedSeqs;
    MR3SeedSites mSeedSites;
    MR3SiteScores mSiteScores;
    mikan::MKRMAWithSites mRNAWithSites;
    MR3SiteFilter mSiteFilter;
    mikan::MKTopNSites mTopNSites;
    MR3RNAScores mRNAScores;

    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_rna_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace mr3as

#endif /* MR3_CORE_HPP_ */
