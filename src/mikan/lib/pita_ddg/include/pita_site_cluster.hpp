#ifndef PITA_SITE_CLUSTER_HPP_
#define PITA_SITE_CLUSTER_HPP_

#include <set>                    // set
#include <map>                    // multimap
#include <utility>                // pair
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"     // MKSiteFilter
#include "pita_score.hpp"         // PITASiteScores
#include "pita_seed_site.hpp"     // PITASeedSites

namespace ptddg {

//
// Identify overlapped sites
//
class PITASiteFilter : public mikan::MKSiteFilter {
public:
    // Define types
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<unsigned, unsigned> TPosPair;

    // Define methods
    PITASiteFilter() : MKSiteFilter() {
        set_gap_len(0);
    }

    void set_gap_len(unsigned pGapLen) { mGapLen = pGapLen; }

private:
    unsigned mGapLen;

    unsigned get_seedtype_precedence(seqan::CharString const &pSeedType);

    void set_intervals(mikan::MKSeedSites &pSeedSites, unsigned pSiteIdx, unsigned &pStartSearch,
                       unsigned &pEndSearch, unsigned &pStartAdd, unsigned &pEndAdd,
                       bool &pSearchOverlap);

};


} // namespace ptddg

#endif /* PITA_SITE_CLUSTER_HPP_ */
