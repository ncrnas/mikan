#ifndef MK_SITE_FILTER_HPP_
#define MK_SITE_FILTER_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_site_score.hpp"        // MKSiteScores
#include "mk_rna_with_sites.hpp"    // MKRMAWithSites

namespace mikan {

//
// Cluster target sites
//
class MKSiteFilter {
public:
    // Define methods
    MKSiteFilter() {
        set_gap_len(0);
    }

    void set_gap_len(unsigned pGapLen) { mGapLen = pGapLen; }

    // Method prototype
    int filter_sites(mikan::MKSeedSites &pSeedSites, mikan::MKRMAWithSites &pRNAWithSites,
                     mikan::MKSiteScores &pSiteScores);

protected:
    // Define method
    virtual float get_precedence(unsigned, mikan::MKSeedSites &, mikan::MKSiteScores &) { return 0; }

    // Method prototypes
    virtual void set_intervals(mikan::MKSeedSites &pSeedSites, unsigned pSiteIdx, unsigned &pStartSearch,
                               unsigned &pEndSearch, unsigned &pStartAdd, unsigned &pEndAdd, bool &pSearchOverlap);

    void mark_overlap(mikan::MKSeedSites &pSeedSites, mikan::TMRNAPosSet &pSortedPos,
                      mikan::MKSiteScores &pSiteScores);

    void sort_sites(mikan::TMRNAPosSet &pSortedPos, mikan::MKSeedSites &pSeedSites,
                    mikan::MKSiteScores &pSiteScores, mikan::TMRNAPosSet &pSortedSeeds);

    // Define types
    typedef std::multimap<float, unsigned> TPosMap;
    typedef std::multimap<float, unsigned>::iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;

    // Define variable
    unsigned mGapLen;

};

} // namespace mikan

#endif /* MK_SITE_FILTER_HPP_ */
