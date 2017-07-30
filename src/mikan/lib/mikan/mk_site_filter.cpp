#include <seqan/seq_io.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"       // MKSiteFilter

using namespace seqan;

namespace mikan {

//
// MKSiteFilter methods
//
int MKSiteFilter::filter_sites(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites &pRNAWithSites,
        mikan::MKSiteScores &pSiteScores) {

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        mark_overlap(pSeedSites, rnaSitePosMap[i], pSiteScores);

    }

    return 0;
}

void MKSiteFilter::mark_overlap(
        mikan::MKSeedSites &pSeedSites,
        mikan::TMRNAPosSet &pSortedPos,
        mikan::MKSiteScores &pSiteScores) {

    IntervalTree<unsigned> tree;
    String<unsigned> results;
    mikan::TMRNAPosSet sortedSeeds;
    unsigned startAdd, endAdd, startSearch, endSearch;
    bool searchOverlap;

    sort_sites(pSortedPos, pSeedSites, pSiteScores, sortedSeeds);

    for (unsigned i = 0; i < length(sortedSeeds); i++) {
        set_intervals(pSeedSites, pSiteScores, sortedSeeds[i], startAdd, endAdd,
                      startSearch, endSearch, searchOverlap);

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

void MKSiteFilter::sort_sites(
        mikan::TMRNAPosSet &pSortedPos,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &pSiteScores,
        mikan::TMRNAPosSet &pSortedSeeds) {

    float preced;
    TPosMap sortedSeeds;

    for (unsigned i = 0; i < length(pSortedPos); i++) {
        preced = get_precedence(pSortedPos[i], pSeedSites, pSiteScores);
        sortedSeeds.insert(TPosPair(preced * (length(pSortedPos) + 1) + i, pSortedPos[i]));
    }

    resize(pSortedSeeds, length(pSortedPos));
    unsigned idx = 0;
    for (TItMap itPos = sortedSeeds.begin(); itPos != sortedSeeds.end(); ++itPos) {
        pSortedSeeds[idx] = (*itPos).second;
        ++idx;
    }

}

void MKSiteFilter::set_intervals(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &,
        unsigned pSiteIdx,
        unsigned &pStartSearch,
        unsigned &pEndSearch,
        unsigned &pStartAdd,
        unsigned &pEndAdd,
        bool &pSearchOverlap) {

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();

    pStartSearch = sitePos[pSiteIdx];
    pEndSearch = pStartSearch + mGapLen;
    pStartAdd = pStartSearch;
    pEndAdd = pEndSearch;
    pSearchOverlap = true;
}

} // namespace mikan
