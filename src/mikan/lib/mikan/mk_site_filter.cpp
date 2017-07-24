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

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();
    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    TPosMap sortedSeeds;

    unsigned pos7, pos1, startPos, endPos;
    CharString seedType;
    bool searchOverlap;

    sort_by_seed_type(pSeedSites, pSortedPos, pCount, sortedSeeds);

    startPos = 0;
    endPos = 0;
    for (TItMap itMap = sortedSeeds.begin(); itMap != sortedSeeds.end(); ++itMap) {
        seedType = seedTypes[(*itMap).second];

        pos7 = sitePos[(*itMap).second] - 1;
        pos1 = sitePos[(*itMap).second] + 7;
        searchOverlap = true;

        if (seedType == "8mer" || seedType == "7mer-A1" || seedType == "7mer-m8" || seedType == "6mer") {
            startPos = pos7;
            endPos = pos1;
            searchOverlap = false;
        } else if (seedType == "GUM" || seedType == "GUT" || seedType == "LP") {
            startPos = pos7 - 1;
            endPos = pos1 + 1;
        } else if (seedType == "BT") {
            startPos = pos7 - 1;
            endPos = pos1 + 2;
        } else if (seedType == "BM") {
            pos7 = sitePos[(*itMap).second];
            pos1 = sitePos[(*itMap).second] + 8;
            startPos = pos7 - 1;
            endPos = pos1;
        }

        clear(results);
        if (searchOverlap) {
            findIntervals(tree, startPos, endPos, results);
        }

        if (length(results) == 0) {
            addInterval(tree, pos7, pos1);
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
    int preced;

    int n = 0;
    for (TItMap itPos = pSortedPos.begin(); itPos != pSortedPos.end(); ++itPos) {
        const CharString &seedType = seedTypes[(*itPos).second];

        if (seedType == "8mer" || seedType == "7mer-A1" || seedType == "7mer-m8") {
            preced = 0;
        } else if (seedType == "6mer") {
            preced = 1;
        } else if (seedType == "GUM") {
            preced = 2;
        } else if (seedType == "GUT") {
            preced = 3;
        } else if (seedType == "BT") {
            preced = 4;
        } else if (seedType == "BM") {
            preced = 5;
        } else if (seedType == "LP") {
            preced = 6;
        } else {
            preced = 7;
        }

        pSortedSeeds.insert(TPosPair(preced * (pCount + 1) + n, (*itPos).second));
        ++n;
    }

}

} // namespace mikan
