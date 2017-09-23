#ifndef TS5_CORE_HPP_
#define TS5_CORE_HPP_

#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_site_filter.hpp"     // MKSiteFilter
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_core_tmpl.hpp"       // MKCoreTmpl
#include "ts5_feature.hpp"        // TS5RawFeatures
#include "ts5_option.hpp"         // TS5Options
#include "ts5_site_score.hpp"     // TS5SiteScores
#include "ts5_seed_site.hpp"      // TS5SeedSites
#include "ts5_rna_score.hpp"      // TS5RNAScores

namespace ts5cs {

//
// TargetScan context score process core
//
typedef mikan::MKCoreTmpl<TS5SeedSeqs, TS5SeedSites, TS5SiteScores, mikan::MKSiteFilter, TS5RNAScores> TS5CoreBase;

class TS5Core : public TS5CoreBase {
public:
    // Define methods
    TS5Core(mikan::MKOptions const &pOpts,
            mikan::TCharSet const &pMiRNAIds,
            mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds,
            mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx,
            mikan::TFinder &pFinder) :
            TS5CoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mClusterSites1 = false;
        mFilterSites = false;
        mClusterSites2 = false;
        mFilterSiteScores = false;
        mSelectTopSites = false;

    }

private:
    virtual int write_site_score(mikan::TCharStr const &pMiRNAId);

    virtual int write_site_score_gff(mikan::TCharStr const &) { return 0; }

    virtual int write_rna_score(mikan::TCharStr const &pMiRNAId);

    virtual int write_rna_score_gff(mikan::TCharStr const &) { return 0; }

    virtual int write_alignment(mikan::TCharStr const &pMiRNAId);

};

} // namespace ts5cs

#endif /* TS5_CORE_HPP_ */
