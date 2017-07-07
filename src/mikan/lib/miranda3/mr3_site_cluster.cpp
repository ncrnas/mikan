#include "mk_typedef.hpp"         // TRNATYPE
#include "mr3_score.hpp"          // MR3MFEScores
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "mr3_site_cluster.hpp"   // MR3Overlap, MR3SortedSitePos
#include <set>                    // set
#include <map>                    // multimap
#include <utility>                // pair
#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;
using namespace mikan;

namespace mr3as {

//
// MR3SiteCluster methods
//
template<class TRNAString>
void MR3SiteCluster<TRNAString>::clear_cluster() {
    mSiteCount = 0;
    mRNAPosSet.clear();
    mSiteMap.clear();
}

template<class TRNAString>
void MR3SiteCluster<TRNAString>::cluster_site_pos(
        MR3SeedSites<TRNAString> &pSeedSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }
        mRNAPosSet.insert(mRNAPos[i]);
        mSiteMap.insert(TPosPair(mRNAPos[i], i));
        ++mSiteCount;
    }
}

//
// MR3Overlap methods
//
template<class TRNAString>
void MR3Overlap<TRNAString>::clear_cluster() {
    mSiteCluster.clear_cluster();
}

template<class TRNAString>
int MR3Overlap<TRNAString>::make_overlapped_pairs(
        MR3SeedSites<TRNAString> &pSeedSites,
        int pGapLen,
        StringSet<String<unsigned> > &pPairs) {
    TItSet itSet;
    TItRetPair ret;
    TItMap itMap;
    TItMap itMap2;

    mSiteCluster.cluster_site_pos(pSeedSites);
    std::set<unsigned> &mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    std::multimap<unsigned, unsigned> startPos;
    TITStartPos itStart;
    int curPos, curIdx, prevPos, prevIdx;

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {
        if (siteMap.count(*itSet) < 2) {
            continue;
        }

        ret = siteMap.equal_range(*itSet);
        for (itMap2 = ret.first; itMap2 != ret.second; ++itMap2) {
            startPos.insert(TPosPair(sitePos[(*itMap2).second], (*itMap2).second));
        }

        prevPos = 0;
        for (itStart = startPos.begin(); itStart != startPos.end(); ++itStart) {
            curPos = (*itStart).first;
            curIdx = (*itStart).second;

            if (prevPos != 0 && ((curPos - prevPos) < (pGapLen + 1))) {
                String<unsigned> tmpPair;
                appendValue(tmpPair, prevIdx);
                appendValue(tmpPair, curIdx);
                appendValue(pPairs, tmpPair);
            }

            prevPos = curPos;
            prevIdx = curIdx;
        }

        startPos.clear();

    }

    return 0;

}

template<class TRNAString>
int MR3Overlap<TRNAString>::filter_overlapped_sites_by_scores(
        MR3SeedSites<TRNAString> &pSeedSites,
        MR3SiteScores<TRNAString> &pSiteScores,
        int pGapLen) {
    StringSet<String<unsigned> > pairs;
    unsigned scorePrev, scoreCur;

    make_overlapped_pairs(pSeedSites, pGapLen, pairs);
    for (unsigned i = 0; i < length(pairs); ++i) {
        if (pSeedSites.mEffectiveSites[pairs[i][0]] == false || pSeedSites.mEffectiveSites[pairs[i][1]] == false) {
            continue;
        }

        scorePrev = (unsigned) pSiteScores.get_align_score(pairs[i][0]);
        scoreCur = (unsigned) pSiteScores.get_align_score(pairs[i][1]);
        if (scorePrev < scoreCur) {
            pSeedSites.mEffectiveSites[pairs[i][0]] = false;
        } else {
            pSeedSites.mEffectiveSites[pairs[i][1]] = false;
        }
    }

    return 0;
}

template<class TRNAString>
int MR3Overlap<TRNAString>::filter_overlapped_sites(MR3SeedSites<TRNAString> &pSeedSites, int pGapLen) {
    StringSet<String<unsigned> > pairs;

    make_overlapped_pairs(pSeedSites, pGapLen, pairs);
    for (unsigned i = 0; length(pairs); ++i) {
        mark_overlapped_sites(pSeedSites, pairs[i][0], pairs[i][1]);
    }

    return 0;
}

template<class TRNAString>
void
MR3Overlap<TRNAString>::mark_overlapped_sites(MR3SeedSites<TRNAString> &pSeedSites, int pPrevIdx, int pCurIdx) {
    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    unsigned precPrev = get_seedtype_precedence(seedTypes[pPrevIdx]);
    unsigned precCur = get_seedtype_precedence(seedTypes[pCurIdx]);

    if (precPrev < precCur) {
        pSeedSites.mEffectiveSites[pCurIdx] = false;
    } else {
        pSeedSites.mEffectiveSites[pPrevIdx] = false;
    }

}


template<class TRNAString>
unsigned MR3Overlap<TRNAString>::get_seedtype_precedence(const CharString &pSeedType) {
    unsigned preced;

    if (pSeedType == "8mer") {
        preced = 0;
    } else if (pSeedType == "7mer") {
        preced = 1;
    } else if (pSeedType == "6mer") {
        preced = 2;
    } else if (pSeedType == "8mer_GUT" || pSeedType == "8mer_GUM") {
        preced = 3;
    } else if (pSeedType == "8mer_GU+") {
        preced = 4;
    } else if (pSeedType == "8mer_MM") {
        preced = 5;
    } else if (pSeedType == "7mer_GUT" || pSeedType == "7mer_GUM") {
        preced = 6;
    } else if (pSeedType == "7mer_GU+") {
        preced = 7;
    } else if (pSeedType == "7mer_MM") {
        preced = 8;
    } else if (pSeedType == "8mer_MMGU") {
        preced = 9;
    } else if (pSeedType == "7mer_MMGU") {
        preced = 10;
    } else if (pSeedType == "8mer_MMGU+") {
        preced = 11;
    } else if (pSeedType == "7mer_MMGU+") {
        preced = 12;
    } else if (pSeedType == "8mer_BT") {
        preced = 13;
    } else if (pSeedType == "7mer_BT") {
        preced = 14;
    } else {
        preced = 15;
    }

    return preced;
}

//
// MR3SortedSitePos methods
//

template<class TRNAString>
void MR3SortedSitePos<TRNAString>::clear_site_pos() {
    clear(mSortedSitePos);
    mSiteCluster.clear_cluster();
}

template<class TRNAString>
int MR3SortedSitePos<TRNAString>::generate_sorted_mrna_pos(
        MR3SeedSites<TRNAString> &pSeedSites) {
    TItMap itMap;
    TItSet itSet;
    TItRetPair ret;
    TITStartPos itStart;
    std::multimap<unsigned, unsigned> startPos;
    int curPos = 0;

    mSiteCluster.cluster_site_pos(pSeedSites);
    std::set<unsigned> &mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();
    int siteCount = mSiteCluster.get_site_count();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();

    resize(mSortedSitePos, siteCount);

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {

        ret = siteMap.equal_range((*itSet));
        for (itMap = ret.first; itMap != ret.second; ++itMap) {
            startPos.insert(TPosPair(sitePos[(*itMap).second], (*itMap).second));
        }

        for (itStart = startPos.begin(); itStart != startPos.end(); ++itStart) {
            mSortedSitePos[curPos] = (*itStart).second;
            ++curPos;
        }

        startPos.clear();
    }

    return 0;

}

// Explicit template instantiation
template
class MR3SiteCluster<TRNATYPE>;

template
class MR3Overlap<TRNATYPE>;

template
class MR3SortedSitePos<TRNATYPE>;

} // namespace mr3as
