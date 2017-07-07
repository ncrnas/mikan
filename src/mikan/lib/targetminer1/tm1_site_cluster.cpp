#include "mk_typedef.hpp"         // TRNATYPE
#include "tm1_seed_site.hpp"      // TM1SeedSites
#include "tm1_site_cluster.hpp"   // TM1Overlap, TM1SortedSitePos
#include <seqan/seq_io.h>

using namespace seqan;
using namespace mikan;

namespace tm1p {

//
// TM1SiteCluster methods
//
template<class TRNAString>
void TM1SiteCluster<TRNAString>::clear_cluster() {
    mSiteCount = 0;
    mRNAPosSet.clear();
    mSiteMap.clear();
}

template<class TRNAString>
void TM1SiteCluster<TRNAString>::cluster_site_pos(
        TM1SeedSites &pSeedSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }
        mRNAPosSet.insert((unsigned) mRNAPos[i]);
        mSiteMap.insert(TPosPair((unsigned) mRNAPos[i], i));
        ++mSiteCount;
    }
}

//
// TM1SortedSitePos methods
//

template<class TRNAString>
void TM1SortedSitePos<TRNAString>::clear_site_pos() {
    clear(mSortedSites);
    clear(mMRNAIDs);

    mSiteCluster.clear_cluster();
}

template<class TRNAString>
int TM1SortedSitePos<TRNAString>::generate_sorted_mrna_pos(
        TM1SeedSites &pSeedSites,
        bool pRemoveOvelaps) {
    TItMMap itMap;
    TItSet itSet;
    TItRetPair ret;
    TITStartPos itStart;
    std::multimap<unsigned, unsigned> startPos;
    int idx1 = 0;
    int idx2 = 0;

    if (pRemoveOvelaps) {
        remove_overlapped_sites(pSeedSites);
    }
    mSiteCluster.cluster_site_pos(pSeedSites);
    std::set<unsigned> &mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();

    resize(mSortedSites, mRNAPosSet.size());
    resize(mMRNAIDs, mRNAPosSet.size());

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {

        mMRNAIDs[idx1] = *itSet;
        ret = siteMap.equal_range((*itSet));
        for (itMap = ret.first; itMap != ret.second; ++itMap) {
            startPos.insert(TPosPair((unsigned) sitePos[(*itMap).second], (*itMap).second));
        }

        resize(mSortedSites[idx1], startPos.size());
        idx2 = 0;
        for (itStart = startPos.begin(); itStart != startPos.end(); ++itStart) {
            mSortedSites[idx1][idx2] = (*itStart).second;
            ++idx2;
        }

        startPos.clear();

        ++idx1;
    }

    return 0;

}

template<class TRNAString>
void TM1SortedSitePos<TRNAString>::remove_overlapped_sites(TM1SeedSites &pSeedSites) {
    TItMMap itMap;
    TItMap itSorted;
    TItSet itSet;
    TItRetPair ret;
    String<unsigned> results;
    std::map<unsigned, unsigned> sortedSites;

    mSiteCluster.cluster_site_pos(pSeedSites);
    std::set<unsigned> &mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();
    unsigned siteId, seedStart, seedEnd;
    bool searchOverlap;

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {
        ret = siteMap.equal_range((*itSet));
        for (itMap = ret.first; itMap != ret.second; ++itMap) {
            sort_by_seed_types(pSeedSites, ret, sortedSites);
        }

        IntervalTree<unsigned> tree;
        searchOverlap = false;
        ret = siteMap.equal_range((*itSet));

        for (itSorted = sortedSites.begin(); itSorted != sortedSites.end(); ++itSorted) {
            siteId = (*itSorted).second;

            seedStart = (unsigned) pSeedSites.get_seed_start_pos(siteId);
            seedEnd = (unsigned) pSeedSites.get_seed_end_pos(siteId) + 1;
            seedEnd = seedStart + 4;

            clear(results);
            if (searchOverlap) {
                findIntervals(tree, seedStart, seedEnd, results);
            }

            if (length(results) == 0) {
                addInterval(tree, seedStart, seedEnd);
            }

            if (searchOverlap && length(results) > 0) {
                pSeedSites.mEffectiveSites[siteId] = false;
            }

            searchOverlap = true;
        }

        sortedSites.clear();

    }

    mSiteCluster.clear_cluster();

}

template<class TRNAString>
void TM1SortedSitePos<TRNAString>::sort_by_pos7(
        TM1SeedSites &pSeedSites,
        TItRetPair &pGroupedSites,
        std::map<unsigned, unsigned> &pSortedSites) {
    TItMMap itMap;
    unsigned siteId, pos7;
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();

    for (itMap = pGroupedSites.first; itMap != pGroupedSites.second; ++itMap) {
        siteId = (*itMap).second;
        pos7 = (unsigned) sitePos[siteId];
        pSortedSites.insert(TPosPair(pos7, siteId));
    }
}

template<class TRNAString>
void TM1SortedSitePos<TRNAString>::sort_by_seed_types(
        TM1SeedSites &pSeedSites,
        TItRetPair &pGroupedSites,
        std::map<unsigned, unsigned> &pSortedSites) {
    TItMMap itMap;
    unsigned siteId, siteOrder, pos7;
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    const StringSet<CharString> &seedTypes = pSeedSites.get_seed_types();
    unsigned maxPos = 0;

    for (itMap = pGroupedSites.first; itMap != pGroupedSites.second; ++itMap) {
        pos7 = (unsigned) sitePos[(*itMap).second];
        if (pos7 > maxPos) {
            maxPos = pos7;
        }
    }

    for (itMap = pGroupedSites.first; itMap != pGroupedSites.second; ++itMap) {
        siteId = (*itMap).second;

        siteOrder = 0;
        if (seedTypes[siteId] == "8mer") {
            siteOrder = 1;
        } else if (seedTypes[siteId] == "7mer-m8") {
            siteOrder = 2;
        } else if (seedTypes[siteId] == "7mer-A1") {
            siteOrder = 3;
        } else if (seedTypes[siteId] == "6mer") {
            siteOrder = 4;
        }

        siteOrder *= (maxPos + 1);
        siteOrder += sitePos[siteId];

        pSortedSites.insert(TPosPair(siteOrder, siteId));
    }
}

// Explicit template instantiation
template
class TM1SiteCluster<TRNATYPE>;

template
class TM1SortedSitePos<TRNATYPE>;

} // namespace tm1p
