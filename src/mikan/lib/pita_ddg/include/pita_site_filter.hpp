#ifndef PITA_SITE_FILTER_HPP_
#define PITA_SITE_FILTER_HPP_

#include <set>                    // set
#include <map>                    // multimap
#include <utility>                // pair
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"     // MKSiteFilter
#include "mk_option.hpp"          // MKOptions
#include "pita_score.hpp"         // PITASiteScores
#include "pita_seed_site.hpp"     // PITASeedSites

namespace ptddg {

//
// Filter overlapped sites
//
class PITASiteFilter : public mikan::MKSiteFilter {
public:
    // Define methods
    PITASiteFilter(mikan::MKOptions const &opts) : MKSiteFilter(opts) {}

private:

    float get_precedence(unsigned pSitePos, mikan::MKSeedSites &pSeedSites,
                         mikan::MKSiteScores &pSiteScores);

};


} // namespace ptddg

#endif /* PITA_SITE_FILTER_HPP_ */
