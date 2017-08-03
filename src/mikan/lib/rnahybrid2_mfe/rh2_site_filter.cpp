#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "rh2_site_score.hpp"    // RH2SiteScores
#include "rh2_site_filter.hpp"   // RH2Overlap, RH2SortedSitePos

using namespace seqan;

namespace rh2mfe {

//
// RH2SiteFilter methods
//
void RH2SiteFilter::init_from_args() {
    mOverlapMethod = mOpts.mOverlapDef;

}

float RH2SiteFilter::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &,
        mikan::MKSiteScores &pSiteScores) {

    float preced = pSiteScores.get_score(pSitePos);

    return preced;
}

void RH2SiteFilter::set_intervals(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &pSiteScores,
        unsigned pSiteIdx,
        unsigned &pStartSearch,
        unsigned &pEndSearch,
        unsigned &pStartAdd,
        unsigned &pEndAdd,
        bool &pSearchOverlap) {

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();

    if (mOverlapMethod == "seed") {
        pStartSearch = sitePos[pSiteIdx];
        pEndSearch = pStartSearch + mikan::SEEDLEN;
    } else {
        pStartSearch = pSiteScores.get_wide_site_start(pSiteIdx);
        pEndSearch = pStartSearch + pSiteScores.get_wide_site_length(pSiteIdx);
    }

    pStartAdd = pStartSearch;
    pEndAdd = pEndSearch;
    pSearchOverlap = true;

}

} // namespace rh2mfe
