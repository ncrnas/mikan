#ifndef MR3_CORE_HPP_
#define MR3_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_core_tmpl.hpp"       // MKCoreTmpl
#include "mr3_option.hpp"         // MR3Options
#include "mr3_site_score.hpp"     // MR3GGDScores
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "mr3_site_filter.hpp"    // MR3SiteFilter
#include "mr3_rna_score.hpp"      // MR3RNAScores

namespace mr3as {

//
// MR3 score process core
//
typedef mikan::MKCoreTmpl<MR3SeedSeqs, MR3SeedSites, MR3SiteScores, MR3SiteFilter, MR3RNAScores> MR3CoreBase;

class MR3Core : public MR3CoreBase {
public:
    // Define methods
    MR3Core(mikan::MKOptions const &pOpts,
            mikan::TCharSet const &pMiRNAIds,
            mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds,
            mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx,
            mikan::TFinder &pFinder) :
            MR3CoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mFindSeedSites = true;
        mFilterSites = false;
        mCalcSiteScore = true;
        mClusterSites1 = false;
        mFilterSiteScores = true;
        mClusterSites2 = true;
        mSelectTopSites = false;
        mClusterSites3 = false;
        mCalcRNAScore = true;
    }

private:
    virtual void prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx);

    virtual void prepare_rna_output(mikan::TCharStr const &pMiRNAId);

    virtual int write_alignment(mikan::TCharStr const &pMiRNAId);

};

} // namespace mr3as

#endif /* MR3_CORE_HPP_ */
