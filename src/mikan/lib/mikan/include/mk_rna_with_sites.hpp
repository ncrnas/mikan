#ifndef MK_RNA_WITH_SITES_HPP_
#define MK_RNA_WITH_SITES_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_site_score.hpp"        // MKSiteScores

namespace mikan {

//
// Cluster target sites
//
class MKRMAWithSites {
public:
    // Define variable
    seqan::String<bool> mEffectiveRNAs;

    // Define methods
    MKRMAWithSites() : mSortVtype("position") {}

    mikan::TMRNAPosSet &get_uniq_mrna_pos_set() { return mUniqRNAPosSet; }

    seqan::StringSet<seqan::String<unsigned> > &get_rna_site_pos_map() { return mRNASitePosMap; }

    // Method prototypes
    void clear_maps();

    void create_mrna_site_map(MKSeedSites &pSeedSites, MKSiteScores &pSiteScores);

    void set_sort_vtype(seqan::CharString &pSortVtype) { mSortVtype = pSortVtype; }

private:
    void create_temp_map(mikan::MKSeedSites &pSeedSites);

    // Define types
    typedef std::set<unsigned> TSet;
    typedef std::pair<float, unsigned> TPosPair;
    typedef std::multimap<float, unsigned> TPosMap;
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<float, unsigned>::iterator TItMap;
    typedef std::pair<TItMap, TItMap> TItMapPair;

    // Define variables
    mikan::TMRNAPosSet mUniqRNAPosSet;
    seqan::StringSet<seqan::String<unsigned> > mRNASitePosMap;
    TSet mUniqSetTemp;
    TPosMap mSiteMapTemp;
    seqan::CharString mSortVtype;

};

} // namespace mikan

#endif /* MK_RNA_WITH_SITES_HPP_ */
