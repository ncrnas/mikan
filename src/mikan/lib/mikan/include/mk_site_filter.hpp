#ifndef MK_SITE_FILTER_HPP_
#define MK_SITE_FILTER_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_rna_with_sites.hpp"    // MKRMAWithSites

namespace mikan {

//
// Cluster target sites
//
class MKSiteFilter {
public:
    // Define methods
    MKSiteFilter() {}

    // Method prototype
    int filter_sites_by_seed_type(mikan::MKSeedSites &pSeedSites, mikan::MKRMAWithSites &pRNAWithSites);

private:
    typedef std::set<unsigned> TSet;
    typedef std::multimap<unsigned, unsigned> TPosMap;
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<TItMap, TItMap> TItMapPair;
    typedef std::pair<unsigned, unsigned> TPosPair;

private:
    virtual void mark_overlap_by_seed_type(mikan::MKSeedSites &pSeedSites, TPosMap &pSortedPos, unsigned pCount);

    virtual void sort_by_seed_type(mikan::MKSeedSites &pSeedSites, TPosMap &pSortedPos, int pCount, TPosMap &sortedSeeds);

};

} // namespace mikan

#endif /* MK_SITE_FILTER_HPP_ */
