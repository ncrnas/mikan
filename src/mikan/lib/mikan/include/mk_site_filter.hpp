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
    // Define methods
    virtual unsigned get_seedtype_precedence(seqan::CharString const &) { return 0; }

    virtual void set_intervals(mikan::MKSeedSites &, unsigned, unsigned &, unsigned &,
                               unsigned &, unsigned &, bool &) {}

    // Method prototypes
    void mark_overlap_by_seed_type(mikan::MKSeedSites &pSeedSites, mikan::TMRNAPosSet &pSortedPos);

    void sort_by_seed_type(mikan::MKSeedSites &pSeedSites, mikan::TMRNAPosSet &pSortedPos,
                           mikan::TMRNAPosSet &pSortedSeeds);

    // Define types
    typedef std::multimap<unsigned, unsigned> TPosMap;
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<unsigned, unsigned> TPosPair;

};

} // namespace mikan

#endif /* MK_SITE_FILTER_HPP_ */
