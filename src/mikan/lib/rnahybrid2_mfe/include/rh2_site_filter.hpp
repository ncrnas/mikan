#ifndef RH2_SITE_FILTER_HPP_
#define RH2_SITE_FILTER_HPP_

#include <set>                   // set
#include <map>                   // multimap
#include <utility>               // pair
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"    // MKSiteFilter
#include "mk_option.hpp"         // MKOptions
#include "rh2_score.hpp"         // RH2SiteScores
#include "rh2_seed_site.hpp"     // RH2SeedSites

namespace rh2mfe {

//
// Filter overlapped sites
//
class RH2SiteFilter : public mikan::MKSiteFilter {
public:
    // Define methods
    RH2SiteFilter(mikan::MKOptions const &opts) : MKSiteFilter(opts) {}

    // Method prototype
    void init_from_args();

private:
    float get_precedence(unsigned pSitePos, mikan::MKSeedSites &pSeedSites,
                         mikan::MKSiteScores &pSiteScores);

    void set_intervals(mikan::MKSeedSites &pSeedSites, mikan::MKSiteScores &pSiteScores, unsigned pSiteIdx,
                       unsigned &pStartSearch, unsigned &pEndSearch, unsigned &pStartAdd, unsigned &pEndAdd,
                       bool &pSearchOverlap);

    seqan::CharString mOverlapMethod;

};

} // namespace rh2mfe

#endif /* RH2_SITE_FILTER_HPP_ */
