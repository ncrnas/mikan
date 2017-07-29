#ifndef MR3_SITE_FILTER_HPP_
#define MR3_SITE_FILTER_HPP_

#include <set>                    // set
#include <map>                    // multimap
#include <utility>                // pair
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"     // MKSiteFilter
#include "mk_option.hpp"          // MKOptions
#include "mr3_score.hpp"          // MR3DDGScores
#include "mr3_seed_site.hpp"      // MR3SeedSites

namespace mr3as {
//
// Filter overlapped sites
//
class MR3SiteFilter : public mikan::MKSiteFilter {
public:
    // Define methods
    MR3SiteFilter(mikan::MKOptions const &opts) : MKSiteFilter(opts) {
            set_gap_len(6);
    }

private:

    float get_precedence(unsigned pSitePos, mikan::MKSeedSites &pSeedSites,
                         mikan::MKSiteScores &pSiteScores);

};

} // namespace mr3as

#endif /* MR3_SITE_FILTER_HPP_ */
