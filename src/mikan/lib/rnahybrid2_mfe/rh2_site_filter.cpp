#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "rh2_score.hpp"         // RH2SiteScores
#include "rh2_site_filter.hpp"  // RH2Overlap, RH2SortedSitePos

using namespace seqan;

namespace rh2mfe {

//
// RH2SiteCluster methods
//
void RH2SiteCluster::clear_cluster() {
    mSiteCount = 0;
    mRNAPosSet.clear();
    mSiteMap.clear();
}

void RH2SiteCluster::cluster_site_pos(
        RH2SeedSites &pSeedSites,
        RH2SiteScores &pScores) {
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
// RH2SiteFilter methods
//
void RH2SiteFilter::init_from_args() {
    mOverlapMethod = mOpts.mOverlapDef;

}

float RH2SiteFilter::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &,
        mikan::MKSiteScores &pSiteScores) {

    float preced = pSiteScores.get_score(pSitePos);

    return preced;
}

void RH2SiteFilter::set_intervals(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &pSiteScores,
        unsigned pSiteIdx,
        unsigned &pStartSearch,
        unsigned &pEndSearch,
        unsigned &pStartAdd,
        unsigned &pEndAdd,
        bool &pSearchOverlap) {

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();

    if (mOverlapMethod == "seed") {
        pStartSearch = sitePos[pSiteIdx];
        pEndSearch = pStartSearch + mikan::SEEDLEN;
    } else {
        pStartSearch = pSiteScores.get_wide_site_start(pSiteIdx);
        pEndSearch = pStartSearch + pSiteScores.get_wide_site_length(pSiteIdx);
    }

    pStartAdd = pStartSearch;
    pEndAdd = pEndSearch;
    pSearchOverlap = true;

}

//
// RH2TopNScore methods
//
void RH2TopNScore::clear_cluster() {
    mTopN = 0;
    mSiteCluster.clear_cluster();
}

int RH2TopNScore::filter_sites(
        RH2SeedSites &pSeedSites,
        RH2SiteScores &pScores,
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

void RH2TopNScore::sort_sites_by_score(RH2SiteScores &pScores, int pPosIdx) {
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

void RH2TopNScore::mark_non_topn_sites(
        RH2SiteScores &pScores,
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

} // namespace rh2mfe
