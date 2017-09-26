#ifndef RH2_CORE_HPP_
#define RH2_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "mk_option.hpp"         // MKOptions
#include "mk_core_tmpl.hpp"      // MKCoreTmpl
#include "rh2_option.hpp"        // RH2Options
#include "rh2_site_score.hpp"    // RH2SiteScores
#include "rh2_seed_site.hpp"     // RH2SeedSites
#include "rh2_site_filter.hpp"   // RH2SiteFilter, RH2TopNSites
#include "rh2_rna_score.hpp"     // RH2RNAScores

namespace rh2mfe {

//
// RNAhybrid MFE score process core
//
typedef mikan::MKCoreTmpl<RH2SeedSeqs, RH2SeedSites, RH2SiteScores, RH2SiteFilter, RH2RNAScores> RH2CoreBase;

class RH2Core : public RH2CoreBase {
public:
    // Define methods
    RH2Core(mikan::MKOptions const &pOpts,
            mikan::TCharSet const &pMiRNAIds,
            mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds,
            mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx,
            mikan::TFinder &pFinder) :
            RH2CoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mClusterSites1 = false;
        mFilterSites = false;
        mClusterSites3 = false;

    }

private:
    virtual void write_site_score_tab(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx);

    virtual void write_site_score_gff(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx);

    virtual void write_rna_score_tab(mikan::TCharStr const &pMiRNAId);

    virtual void write_rna_score_gff(mikan::TCharStr const &) { return; }

    virtual int write_alignment(mikan::TCharStr const &pMiRNAId);

};

} // namespace rh2mfe

#endif /* RH2_CORE_HPP_ */
