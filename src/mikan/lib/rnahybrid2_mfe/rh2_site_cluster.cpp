#include "mk_typedef.hpp"        // TRNATYPE
#include "rh2_score.hpp"         // RH2MFEScores
#include "rh2_site_cluster.hpp"  // RH2Overlap, RH2SortedSitePos

using namespace seqan;
using namespace mikan;

namespace rh2mfe {

//
// RH2SiteCluster methods
//
template<class TRNAString>
void RH2SiteCluster<TRNAString>::clear_cluster() {
    mSiteCount = 0;
    mRNAPosSet.clear();
    mSiteMap.clear();
}

template<class TRNAString>
void RH2SiteCluster<TRNAString>::cluster_site_pos(
        RH2SeedSites<TRNAString> &pSeedSites,
        RH2MFEScores<TRNAString> &pScores) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pScores.mEffectiveSites[i]) {
            continue;
        }
        mRNAPosSet.insert((unsigned) mRNAPos[i]);
        mSiteMap.insert(TPosPair((unsigned) mRNAPos[i], i));
        ++mSiteCount;
    }
}

//
// RH2Overlap methods
//
template<class TRNAString>
void RH2Overlap<TRNAString>::clear_cluster() {
    mSiteCluster.clear_cluster();
}

template<class TRNAString>
int RH2Overlap<TRNAString>::filter_overlapped_sites(
        RH2SeedSites<TRNAString> &pSeedSites,
        RH2MFEScores<TRNAString> &pScores,
        CharString &pOverlapDef) {
    TItSet itSet;

    mSiteCluster.cluster_site_pos(pSeedSites, pScores);
    std::set<unsigned> &mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {
        if (siteMap.count(*itSet) > 1) {
            find_overlapped_sites(pSeedSites, pScores, *itSet, pOverlapDef);
        }
    }

    return 0;

}

template<class TRNAString>
void RH2Overlap<TRNAString>::find_overlapped_sites(
        RH2SeedSites<TRNAString> &pSeedSites,
        RH2MFEScores<TRNAString> &pScores,
        int pPosIdx,
        CharString &pOverlapDef) {
    TItMap itMap;
    TItRetPair ret;
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();
    std::multimap<unsigned, unsigned> startPos;
    const seqan::String<unsigned> &sitePos = pSeedSites.get_site_pos();
    unsigned stpos;

    ret = siteMap.equal_range((unsigned) pPosIdx);
    for (itMap = ret.first; itMap != ret.second; ++itMap) {
        if (pOverlapDef == "seed") {
            stpos = sitePos[(*itMap).second];
        } else {
            stpos = (unsigned) pScores.get_hit_start((*itMap).second);
        }
        startPos.insert(TPosPair(stpos, (*itMap).second));
    }

    cluster_overlapped_sites(pSeedSites, pScores, startPos, pOverlapDef);
}

template<class TRNAString>
void RH2Overlap<TRNAString>::cluster_overlapped_sites(
        RH2SeedSites<TRNAString> &pSeedSites,
        RH2MFEScores<TRNAString> &pScores,
        std::multimap<unsigned, unsigned> &pStartPos,
        CharString &pOverlapDef) {
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    TITStartPos itStart;
    std::set<unsigned> olCluster;
    unsigned prevStartPos = 0;
    unsigned prevMRNAPos = 0;

    for (itStart = pStartPos.begin(); itStart != pStartPos.end(); ++itStart) {
        if (sitePos[(*itStart).second] < prevStartPos) {
            olCluster.insert(prevMRNAPos);
            olCluster.insert((*itStart).second);
        } else if (olCluster.size() > 0) {
            mark_overlapped_sites(pSeedSites, pScores, olCluster, pOverlapDef);;
            olCluster.clear();
        }
        if (pOverlapDef == "seed") {
            prevStartPos = (*itStart).first + RH2SeedSites<TRNAString>::SEED_LEN;
        } else {
            prevStartPos = (*itStart).first + pScores.get_hit_length((*itStart).second);
        }
        prevMRNAPos = (*itStart).second;
    }
    if (olCluster.size() > 0) {
        mark_overlapped_sites(pSeedSites, pScores, olCluster, pOverlapDef);

        olCluster.clear();
    }

}

template<class TRNAString>
void RH2Overlap<TRNAString>::mark_overlapped_sites(
        RH2SeedSites<TRNAString> &pSeedSites,
        RH2MFEScores<TRNAString> &pScores,
        std::set<unsigned> &pOlCluster,
        CharString &pOverlapDef) {
    TItSet itSet;
    unsigned mPos;
    std::multimap<unsigned, unsigned> startPos;
    const seqan::String<unsigned> &sitePos = pSeedSites.get_site_pos();
    unsigned stpos1;
    unsigned stpos2;

    mPos = get_pos_with_best_mfe(pScores, pOlCluster);

    if (pOlCluster.size() > 2) {
        for (itSet = pOlCluster.begin(); itSet != pOlCluster.end(); ++itSet) {
            if (mPos == (*itSet) || !pScores.mEffectiveSites[(*itSet)]) {
                continue;
            }
            if (pOverlapDef == "seed") {
                stpos1 = sitePos[mPos];
                stpos2 = sitePos[(*itSet)];
            } else {
                stpos1 = (unsigned) pScores.get_hit_start(mPos);
                stpos2 = (unsigned) pScores.get_hit_start((*itSet));
            }
            startPos.insert(TPosPair(stpos1, (mPos)));
            startPos.insert(TPosPair(stpos2, (*itSet)));
            cluster_overlapped_sites(pSeedSites, pScores, startPos, pOverlapDef);
            startPos.clear();
        }

        for (itSet = pOlCluster.begin(); itSet != pOlCluster.end(); ++itSet) {
            if (mPos == (*itSet) || !pScores.mEffectiveSites[(*itSet)]) {
                continue;
            }
            if (pOverlapDef == "seed") {
                stpos2 = sitePos[(*itSet)];
            } else {
                stpos2 = (unsigned) pScores.get_hit_start((*itSet));
            }
            startPos.insert(TPosPair(stpos2, (*itSet)));
        }
        cluster_overlapped_sites(pSeedSites, pScores, startPos, pOverlapDef);
        startPos.clear();
    } else {
        for (itSet = pOlCluster.begin(); itSet != pOlCluster.end(); ++itSet) {
            if (mPos == (*itSet)) {
                continue;
            }
            pScores.mEffectiveSites[*itSet] = false;
        }
    }

}

template<class TRNAString>
unsigned RH2Overlap<TRNAString>::get_pos_with_best_mfe(
        RH2MFEScores<TRNAString> &pScores,
        std::set<unsigned> &pOlCluster) {
    TItSet itSet;
    unsigned mPos = 0;
    float bestMfe = 65000;

    for (itSet = pOlCluster.begin(); itSet != pOlCluster.end(); ++itSet) {
        if (bestMfe > pScores.get_score(*itSet)) {
            bestMfe = pScores.get_score(*itSet);
            mPos = *itSet;
        }
    }

    return mPos;

}

//
// RH2TopNScore methods
//
template<class TRNAString>
void RH2TopNScore<TRNAString>::clear_cluster() {
    mTopN = 0;
    mSiteCluster.clear_cluster();
}

template<class TRNAString>
int RH2TopNScore<TRNAString>::filter_sites(
        RH2SeedSites<TRNAString> &pSeedSites,
        RH2MFEScores<TRNAString> &pScores,
        int pMaxHits) {
    TItSet itSet;

    if (pMaxHits == 0) {
        return 0;
    }
    mTopN = pMaxHits;

    mSiteCluster.cluster_site_pos(pSeedSites, pScores);
    std::set<unsigned> &mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {
        if (siteMap.count(*itSet) > 1) {
            sort_sites_by_score(pScores, *itSet);
        }
    }

    return 0;

}

template<class TRNAString>
void RH2TopNScore<TRNAString>::sort_sites_by_score(RH2MFEScores<TRNAString> &pScores, int pPosIdx) {
    TItMap itMap;
    TItRetPair ret;
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();
    std::multimap<float, unsigned> sortedSites;

    ret = siteMap.equal_range((unsigned) pPosIdx);
    for (itMap = ret.first; itMap != ret.second; ++itMap) {
        sortedSites.insert(TScorePair(pScores.get_score((*itMap).second), (*itMap).second));
    }

    mark_non_topn_sites(pScores, sortedSites);
}

template<class TRNAString>
void RH2TopNScore<TRNAString>::mark_non_topn_sites(
        RH2MFEScores<TRNAString> &pScores,
        std::multimap<float, unsigned> &pSortedSites) {
    TITStartScore itStart;
    int siteCount = 0;

    for (itStart = pSortedSites.begin(); itStart != pSortedSites.end(); ++itStart) {
        if (siteCount < mTopN) {
            ++siteCount;
            continue;
        }

        pScores.mEffectiveSites[(*itStart).second] = false;
    }

}

//
// RH2SortedSitePos methods
//

template<class TRNAString>
void RH2SortedSitePos<TRNAString>::clear_site_pos() {
    clear(mSortedSitePos);
    mSiteCluster.clear_cluster();
}

template<class TRNAString>
int RH2SortedSitePos<TRNAString>::generate_sorted_mrna_pos(
        RH2SeedSites<TRNAString> &pSeedSites,
        RH2MFEScores<TRNAString> &pScores) {
    TItMap itMap;
    TItSet itSet;
    TItRetPair ret;
    TITStartPos itStart;
    std::multimap<unsigned, unsigned> startPos;
    int curPos = 0;

    mSiteCluster.cluster_site_pos(pSeedSites, pScores);
    std::set<unsigned> &mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned> &siteMap = mSiteCluster.get_mrna_pos_map();
    int siteCount = mSiteCluster.get_site_count();

    resize(mSortedSitePos, siteCount);

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {

        ret = siteMap.equal_range((*itSet));
        for (itMap = ret.first; itMap != ret.second; ++itMap) {
            startPos.insert(TPosPair((unsigned) pScores.get_hit_start((*itMap).second), (*itMap).second));
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
class RH2SiteCluster<TRNATYPE>;

template
class RH2Overlap<TRNATYPE>;

template
class RH2TopNScore<TRNATYPE>;

template
class RH2SortedSitePos<TRNATYPE>;

} // namespace rh2mfe
