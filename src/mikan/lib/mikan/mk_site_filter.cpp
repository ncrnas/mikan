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

    if (!mUseFilter) {
        return 0;
    }

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();

    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i] || (length(rnaSitePosMap[i]) <= 1)) {
            continue;
        }

        mark_sites(pSeedSites, rnaSitePosMap[i], pSiteScores);

    }

    return 0;
}

void MKSiteFilter::mark_sites(
        mikan::MKSeedSites &pSeedSites,
        mikan::TMRNAPosSet &pSortedPos,
        mikan::MKSiteScores &pSiteScores) {

    IntervalTree<unsigned> tree;
    String<unsigned> results;
    mikan::TMRNAPosSet sortedSites;
    unsigned startAdd, endAdd, startSearch, endSearch;
    bool searchOverlap;

    sort_sites(pSortedPos, pSeedSites, pSiteScores, sortedSites);

    for (unsigned i = 0; i < length(sortedSites); i++) {
        set_intervals(pSeedSites, pSiteScores, sortedSites[i], startAdd, endAdd,
                      startSearch, endSearch, searchOverlap);

        clear(results);
        if (searchOverlap) {
            findIntervals(tree, startSearch, endSearch, results);
        }

        if (length(results) == 0) {
            addInterval(tree, startAdd, endAdd);
        }

        if (searchOverlap && length(results) > 0) {
            pSeedSites.mEffectiveSites[sortedSites[i]] = false;
        }
    }

}

void MKSiteFilter::sort_sites(
        mikan::TMRNAPosSet &pSortedPos,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &pSiteScores,
        mikan::TMRNAPosSet &pSortedSites) {

    float preced;
    TPosMap sortedSeeds;

    for (unsigned i = 0; i < length(pSortedPos); i++) {
        preced = get_precedence(pSortedPos[i], pSeedSites, pSiteScores);
        sortedSeeds.insert(TPosPair(preced * (length(pSortedPos) + 1) + i, pSortedPos[i]));
    }

    resize(pSortedSites, length(pSortedPos));
    unsigned idx = 0;
    for (TItMap itPos = sortedSeeds.begin(); itPos != sortedSeeds.end(); ++itPos) {
        pSortedSites[idx] = (*itPos).second;
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
    pEndSearch = pStartSearch + mOverlapLen;
    pStartAdd = pStartSearch;
    pEndAdd = pEndSearch;
    pSearchOverlap = true;
}

//
// MKTopNSites methods
//

void MKTopNSites::init_from_args() {
    mTopN = mOpts.mMaxHits;
    if (mTopN == 0) {
        mUseFilter = false;
    }

}

void MKTopNSites::mark_sites(
        mikan::MKSeedSites &pSeedSites,
        mikan::TMRNAPosSet &pSortedPos,
        mikan::MKSiteScores &pSiteScores) {

    int siteCount = 0;
    mikan::TMRNAPosSet sortedSites;

    sort_sites(pSortedPos, pSeedSites, pSiteScores, sortedSites);

    for (unsigned i = 0; i < length(sortedSites); i++) {
        if (siteCount < mTopN) {
            ++siteCount;
            continue;
        }

        pSeedSites.mEffectiveSites[sortedSites[i]] = false;
        pSiteScores.mEffectiveSites[sortedSites[i]] = false;
    }

}

float MKTopNSites::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &,
        mikan::MKSiteScores &pSiteScores) {

    float preced = pSiteScores.get_score(pSitePos);

    return preced;
}

} // namespace mikan
