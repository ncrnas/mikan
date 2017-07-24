#include <seqan/seq_io.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"       // MKSiteFilter

using namespace seqan;

namespace mikan {

//
// MKSiteFilter methods
//
int MKSiteFilter::filter_sites_by_seed_type(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites &pRNAWithSites) {

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();
    TSet uniqRNASet = pRNAWithSites.get_uniq_mrna_set();
    TPosMap rnaSiteMap = pRNAWithSites.get_rna_site_map();

    for (TItSet itSet = uniqRNASet.begin(); itSet != uniqRNASet.end(); ++itSet) {
        unsigned count = rnaSiteMap.count(*itSet);
        if (count < 2) {
            continue;
        }

        TItMapPair itPair = rnaSiteMap.equal_range(*itSet);
        TPosMap sortedPos;
        for (TItMap itMap = itPair.first; itMap != itPair.second; ++itMap) {
            sortedPos.insert(TPosPair(static_cast<unsigned>(sitePos[(*itMap).second]),
                                      (*itMap).second));
        }

        mark_overlap_by_seed_type(pSeedSites, sortedPos, count);
    }

    return 0;
}

void MKSiteFilter::mark_overlap_by_seed_type(
        mikan::MKSeedSites &pSeedSites,
        TPosMap &pSortedPos,
        unsigned pCount) {

    IntervalTree<unsigned> tree;
    String<unsigned> results;
    TPosMap sortedSeeds;
    unsigned startAdd, endAdd, startSearch, endSearch;
    bool searchOverlap;

    sort_by_seed_type(pSeedSites, pSortedPos, pCount, sortedSeeds);

    for (TItMap itMap = sortedSeeds.begin(); itMap != sortedSeeds.end(); ++itMap) {

        set_intervals(pSeedSites, (*itMap).second, startAdd, endAdd, startSearch, endSearch, searchOverlap);

        clear(results);
        if (searchOverlap) {
            findIntervals(tree, startSearch, endSearch, results);
        }

        if (length(results) == 0) {
            addInterval(tree, startAdd, endAdd);
        }

        if (searchOverlap && length(results) > 0) {
            pSeedSites.mEffectiveSites[(*itMap).second] = false;
        }
    }

}

void MKSiteFilter::sort_by_seed_type(
        mikan::MKSeedSites &pSeedSites,
        TPosMap &pSortedPos,
        int pCount,
        TPosMap &pSortedSeeds) {

    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    unsigned preced;

    int n = 0;
    for (TItMap itPos = pSortedPos.begin(); itPos != pSortedPos.end(); ++itPos) {
        preced = get_seedtype_precedence(seedTypes[(*itPos).second]);
        pSortedSeeds.insert(TPosPair(preced * (pCount + 1) + n, (*itPos).second));
        ++n;
    }

}

} // namespace mikan
