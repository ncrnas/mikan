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

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        mark_overlap_by_seed_type(pSeedSites, rnaSitePosMap[i]);

    }

    return 0;
}

void MKSiteFilter::mark_overlap_by_seed_type(
        mikan::MKSeedSites &pSeedSites,
        mikan::TMRNAPosSet &pSortedPos) {

    IntervalTree<unsigned> tree;
    String<unsigned> results;
    mikan::TMRNAPosSet sortedSeeds;
    unsigned startAdd, endAdd, startSearch, endSearch;
    bool searchOverlap;

    sort_by_seed_type(pSeedSites, pSortedPos, sortedSeeds);

    for (unsigned i = 0; i < length(sortedSeeds); i++) {
        set_intervals(pSeedSites, sortedSeeds[i], startAdd, endAdd, startSearch, endSearch, searchOverlap);

        clear(results);
        if (searchOverlap) {
            findIntervals(tree, startSearch, endSearch, results);
        }

        if (length(results) == 0) {
            addInterval(tree, startAdd, endAdd);
        }

        if (searchOverlap && length(results) > 0) {
            pSeedSites.mEffectiveSites[sortedSeeds[i]] = false;
        }
    }

}

void MKSiteFilter::sort_by_seed_type(
        mikan::MKSeedSites &pSeedSites,
        mikan::TMRNAPosSet &pSortedPos,
        mikan::TMRNAPosSet &pSortedSeeds) {

    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    unsigned preced;
    TPosMap sortedSeeds;

    for (unsigned i = 0; i < length(pSortedPos); i++) {
        preced = get_seedtype_precedence(seedTypes[pSortedPos[i]]);
        sortedSeeds.insert(TPosPair(preced * (length(pSortedPos) + 1) + i, pSortedPos[i]));
    }

    resize(pSortedSeeds, length(pSortedPos));
    unsigned idx = 0;
    for (TItMap itPos = sortedSeeds.begin(); itPos != sortedSeeds.end(); ++itPos) {
        pSortedSeeds[idx] = (*itPos).second;
        ++idx;
    }

}

} // namespace mikan
