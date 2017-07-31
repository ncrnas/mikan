#ifndef TM1_SITE_FILTER_HPP_
#define TM1_SITE_FILTER_HPP_

#include <set>                    // set
#include <map>                    // multimap
#include <utility>                // pair
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"     // MKSiteFilter
#include "tm1_seed_site.hpp"      // TM1SeedSites

namespace tm1p {

//
// Filter overlapped sites
//
class TM1SiteFilter : public mikan::MKSiteFilter {
public:
    // Define methods
    TM1SiteFilter(mikan::MKOptions const &opts) : MKSiteFilter(opts) {
        set_overlap_len(4);
    }

private:
    float get_precedence(unsigned pSitePos, mikan::MKSeedSites &pSeedSites,
                         mikan::MKSiteScores &pSiteScores);

};

} // namespace tm1p

#endif /* TM1_SITE_FILTER_HPP_ */
