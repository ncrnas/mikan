#ifndef MK_RNA_WITH_SITES_HPP_
#define MK_RNA_WITH_SITES_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_site.hpp"         // MKSeedSites

namespace mikan {

//
// Cluster target sites
//
class MKRMAWithSites {
public:
    // Define types
    typedef std::set<unsigned> TSet;
    typedef std::pair<unsigned, unsigned> TPosPair;
    typedef std::multimap<unsigned, unsigned> TPosMap;
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<TItMap, TItMap> TItMapPair;

    // Define methods
    MKRMAWithSites() {
        mEffectiveSiteCount = 0;
    }

    TSet &get_uniq_mrna_set() { return mUniqRNASet; }

    TPosMap &get_rna_site_map() { return mRNASiteMap; }

    unsigned get_effective_site_count() { return mEffectiveSiteCount; }

    seqan::StringSet<seqan::String<unsigned> > &get_sorted_mrna_pos() { return mSortedMRNAPos; }

    // Method prototypes
    void clear_sets();

    void cluster_sites(MKSeedSites &pSeedSites);

    void sort_mrna_pos(MKSeedSites &pSeedSites);

private:

    // Define variables
    unsigned mEffectiveSiteCount;
    std::set<unsigned> mUniqRNASet;
    std::multimap<unsigned, unsigned> mRNASiteMap;

    seqan::StringSet<seqan::String<unsigned> > mSortedMRNAPos;
};

} // namespace mikan

#endif /* MK_RNA_WITH_SITES_HPP_ */
