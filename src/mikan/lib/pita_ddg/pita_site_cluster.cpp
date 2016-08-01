#include <mikan/lib/pita_ddg/include/pita_inst_template.hpp> // TRNATYPE
#include <mikan/lib/pita_ddg/include/pita_score.hpp>         // PITAMFEScores
#include <mikan/lib/pita_ddg/include/pita_site_cluster.hpp>  // PITAOverlap, PITASortedSitePos

using namespace seqan;

namespace ptddg{

//
// PITASiteCluster methods
//
template <class TRNAString>
void PITASiteCluster<TRNAString>::clear_cluster()
{
    mSiteCount = 0;
    mRNAPosSet.clear();
    mSiteMap.clear();
}

template <class TRNAString>
void PITASiteCluster<TRNAString>::cluster_site_pos(
        PITASeedSites<TRNAString> &pSeedSites)
{
    const String<unsigned>& mRNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(mRNAPos); ++i)
    {
        if (!pSeedSites.mEffectiveSites[i])
        {
            continue;
        }
        mRNAPosSet.insert((unsigned)mRNAPos[i]);
        mSiteMap.insert(TPosPair((unsigned)mRNAPos[i], i));
        ++mSiteCount;
    }
}

//
// PITAOverlap methods
//
template <class TRNAString>
void PITAOverlap<TRNAString>::clear_cluster()
{
    mSiteCluster.clear_cluster();
}

template <class TRNAString>
int PITAOverlap<TRNAString>::filter_overlapped_sites(PITASeedSites<TRNAString> &pSeedSites, int pGapLen)
{
    TItSet itSet;
    TItRetPair ret;
    TItMap itMap2;

    mSiteCluster.cluster_site_pos(pSeedSites);
    std::set<unsigned>& mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned>& siteMap = mSiteCluster.get_mrna_pos_map();
    const String<unsigned>& sitePos = pSeedSites.get_site_pos();
    std::multimap<unsigned, unsigned> startPos;
    TITStartPos itStart;
    int curPos, curIdx, prevPos, prevIdx;

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet)
    {
        if (siteMap.count(*itSet) < 2)
        {
            continue;
        }

        ret = siteMap.equal_range(*itSet);
        for (itMap2 = ret.first; itMap2 != ret.second; ++itMap2)
        {
            startPos.insert(TPosPair((unsigned)sitePos[(*itMap2).second], (*itMap2).second));
        }

        prevPos = 0;
        prevIdx = 0;
        for (itStart = startPos.begin(); itStart != startPos.end(); ++itStart)
        {
            curPos = (*itStart).first;
            curIdx = (*itStart).second;

            if (prevPos != 0 && ((curPos - prevPos) < (pGapLen + 1)))
            {
                mark_overlapped_sites(pSeedSites, prevIdx, curIdx);
            }

            prevPos = curPos;
            prevIdx = curIdx;
        }

        startPos.clear();

    }

    return 0;

}

template <class TRNAString>
void PITAOverlap<TRNAString>::mark_overlapped_sites(PITASeedSites<TRNAString> &pSeedSites, int pPrevIdx, int pCurIdx)
{
    StringSet<CharString> const& seedTypes = pSeedSites.get_seed_types();
    unsigned precPrev = get_seedtype_precedence(seedTypes[pPrevIdx]);
    unsigned precCur = get_seedtype_precedence(seedTypes[pCurIdx]);

    if (precPrev < precCur)
    {
        pSeedSites.mEffectiveSites[pCurIdx] = false;
    }
    else
    {
        pSeedSites.mEffectiveSites[pPrevIdx] = false;
    }

}


template <class TRNAString>
unsigned PITAOverlap<TRNAString>::get_seedtype_precedence(const CharString &pSeedType)
{
    unsigned preced;

    if (pSeedType == "8mer")
    {
        preced = 0;
    }
    else if (pSeedType == "7mer")
    {
        preced = 1;
    }
    else if (pSeedType == "6mer")
    {
        preced = 2;
    }
    else if (pSeedType == "8mer_GUT" || pSeedType == "8mer_GUM")
    {
        preced = 3;
    }
    else if (pSeedType == "8mer_GU+")
    {
        preced = 4;
    }
    else if (pSeedType == "8mer_MM")
    {
        preced = 5;
    }
    else if (pSeedType == "7mer_GUT" || pSeedType == "7mer_GUM")
    {
        preced = 6;
    }
    else if (pSeedType == "7mer_GU+")
    {
        preced = 7;
    }
    else if (pSeedType == "7mer_MM")
    {
        preced = 8;
    }
    else if (pSeedType == "8mer_MMGU")
    {
        preced = 9;
    }
    else if (pSeedType == "7mer_MMGU")
    {
        preced = 10;
    }
    else if (pSeedType == "8mer_MMGU+")
    {
        preced = 11;
    }
    else if (pSeedType == "7mer_MMGU+")
    {
        preced = 12;
    }
    else
    {
        preced = 13;
    }

    return preced;
}

//
// PITASortedSitePos methods
//

template <class TRNAString>
void PITASortedSitePos<TRNAString>::clear_site_pos()
{
    clear(mSortedSitePos);
    mSiteCluster.clear_cluster();
}

template <class TRNAString>
int PITASortedSitePos<TRNAString>::generate_sorted_mrna_pos(
        PITASeedSites<TRNAString> &pSeedSites)
{
    TItMap itMap;
    TItSet itSet;
    TItRetPair ret;
    TITStartPos itStart;
    std::multimap<unsigned, unsigned> startPos;
    int curPos = 0;

    mSiteCluster.cluster_site_pos(pSeedSites);
    std::set<unsigned>& mRNAPosSet = mSiteCluster.get_mrna_pos_set();
    std::multimap<unsigned, unsigned>& siteMap = mSiteCluster.get_mrna_pos_map();
    int siteCount = mSiteCluster.get_site_count();
    const String<unsigned>& sitePos = pSeedSites.get_site_pos();

    resize(mSortedSitePos, siteCount);

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet)
    {

        ret = siteMap.equal_range((*itSet));
        for (itMap = ret.first; itMap != ret.second; ++itMap)
        {
            startPos.insert(TPosPair((unsigned)sitePos[(*itMap).second], (*itMap).second));
        }

        for (itStart = startPos.begin(); itStart != startPos.end(); ++itStart)
        {
            mSortedSitePos[curPos] = (*itStart).second;
            ++curPos;
        }

        startPos.clear();
    }

    return 0;

}

// Explicit template instantiation
template class PITASiteCluster<TRNATYPE>;
template class PITAOverlap<TRNATYPE>;
template class PITASortedSitePos<TRNATYPE>;

} // namespace ptddg
