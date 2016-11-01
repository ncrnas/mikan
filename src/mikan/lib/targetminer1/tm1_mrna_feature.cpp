#include <tm1_inst_template.hpp> // TRNATYPE
#include <tm1_mrna_feature.hpp>  // TM1ScaledFeatures, TM1RawMRNAures, TM1MRNASeedType, TM1MRNASiteCount,

using namespace seqan;

namespace tm1p{

//
// TM1ScaledFeatures methods
//
template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::init_max_vals()
{

    resize(mMaxFeatVals, FEATURE_NUM);
    mMaxFeatVals[0] = 22.0f;
    mMaxFeatVals[1] = 16.0f;
    mMaxFeatVals[2] = 4.0f;
    mMaxFeatVals[3] = 4.0f;
    mMaxFeatVals[4] = 4.0f;
    mMaxFeatVals[5] = 5.521905f;
    mMaxFeatVals[6] = 5.098095f;
    mMaxFeatVals[7] = 4.668748f;
    mMaxFeatVals[8] = 6.165895f;
    mMaxFeatVals[9] = 3.253333f;
    mMaxFeatVals[10] = 1.092381f;
    mMaxFeatVals[11] = 1.194862f;
    mMaxFeatVals[12] = 1.56094f;
    mMaxFeatVals[13] = 1.265512f;
    mMaxFeatVals[14] = 0.875394f;
    mMaxFeatVals[15] = 0.95453f;
    mMaxFeatVals[16] = 1.444317f;
    mMaxFeatVals[17] = 2.26501f;
    mMaxFeatVals[18] = 1.229357f;
    mMaxFeatVals[19] = 1.152942f;
    mMaxFeatVals[20] = 1.446197f;
    mMaxFeatVals[21] = 0.671929f;
    mMaxFeatVals[22] = 0.9738289999999999f;
    mMaxFeatVals[23] = 1.0f;
    mMaxFeatVals[24] = 1.0f;
    mMaxFeatVals[25] = 4.0f;
    mMaxFeatVals[26] = 2.0f;
    mMaxFeatVals[27] = 1.0f;
    mMaxFeatVals[28] = 2.0f;
    mMaxFeatVals[29] = 1.0f;
}

template <class TRNAString>
float TM1ScaledFeatures<TRNAString>::scale_feat_val(unsigned pIdx, float pVal)
{
    float val;

    if(pVal == 0)
    {
        val = mLower;
    }
    else if(pVal == mMaxFeatVals[pIdx])
    {
        val = mUpper;
    }
    else
    {
        val = mLower + (mUpper - mLower) * (pVal / mMaxFeatVals[pIdx]);
    }

    return val;
}


template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::clear_features()
{
    clear(mMRNAIDs);
    clear(mScaledFeatures);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::resize_features(unsigned pSize)
{
    resize(mMRNAIDs, pSize);
    resize(mScaledFeatures, pSize);
    for (unsigned i = 0; i < pSize; ++i)
    {
        resize(mScaledFeatures[i], FEATURE_NUM);
    }
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::print_features()
{
    for (unsigned i = 0; i < length(mMRNAIDs); ++i)
    {
        for (unsigned j = 0; j < FEATURE_NUM; ++j)
        {
            std::cout << j << ":";
            std::cout << mScaledFeatures[i][j];
            if (j < FEATURE_NUM - 1)
            {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }
}

template <class TRNAString>
int TM1ScaledFeatures<TRNAString>::scale_features(
        TM1SortedSitePos<TRNAString>& pSortedSites,
        TM1MRNASeedType<TRNAString>& pSeedTypes,
        TM1MRNAAURich<TRNAString>& pAURich,
        TM1MRNASingleFreq<TRNAString>& pSingleFreqs,
        TM1MRNASingleFreqFlank<TRNAString>& pSingleFreqFlanks,
        TM1MRNADiFreq<TRNAString>& pDiFreqs,
        TM1MRNADiFreqFlank<TRNAString>& pDiFreqFlanks,
        TM1MRNASingleMatch<TRNAString>& pSingleMatches,
        TM1MRNATwoConsecMatch<TRNAString>& pTwoMatches)
{
    const String<unsigned>& mRNAIDs = pSortedSites.get_mrna_ids();

    resize_features((unsigned)length(mRNAIDs));
    for (unsigned i = 0; i < length(mRNAIDs); ++i)
    {
        mMRNAIDs[i] = mRNAIDs[i];
    }

    for (unsigned i = 0; i < length(mMRNAIDs); ++i)
    {
        scale_seed_type(i, pSeedTypes);
        scale_au_rich(i, pAURich);
        scale_single_freq(i, pSingleFreqs);
        scale_freq_flank(i, pSingleFreqFlanks);
        scale_di_freq(i, pDiFreqs);
        scale_di_freq_flank(i, pDiFreqFlanks);
        scale_single_match(i, pSingleMatches);
        scale_two_match(i, pTwoMatches);
    }

    return 0;
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_seed_type(unsigned pIdx, TM1MRNASeedType<TRNAString>& pSeedTypes)
{
    mScaledFeatures[pIdx][0] = scale_feat_val(0, pSeedTypes.mNum6mer[pIdx]);
    mScaledFeatures[pIdx][3] = scale_feat_val(3, pSeedTypes.mNum8mer[pIdx]);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_au_rich(unsigned pIdx, TM1MRNAAURich<TRNAString>& pAURich)
{
    mScaledFeatures[pIdx][1] = scale_feat_val(1, pAURich.mAU6mer[pIdx]);
    mScaledFeatures[pIdx][2] = scale_feat_val(2, pAURich.mAU7merm8[pIdx]);
    mScaledFeatures[pIdx][4] = scale_feat_val(4, pAURich.mAU8mer[pIdx]);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_single_freq(unsigned pIdx, TM1MRNASingleFreq<TRNAString>& pSingleFreqs)
{
    mScaledFeatures[pIdx][5] = scale_feat_val(5, pSingleFreqs.mSeedFreqU[pIdx]);
    mScaledFeatures[pIdx][6] = scale_feat_val(6, pSingleFreqs.mSeedFreqC[pIdx]);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_freq_flank(
        unsigned pIdx,
        TM1MRNASingleFreqFlank<TRNAString>& pSingleFreqFlanks)
{
    mScaledFeatures[pIdx][7] = scale_feat_val(7, pSingleFreqFlanks.m3PFreqA[pIdx]);
    mScaledFeatures[pIdx][8] = scale_feat_val(8, pSingleFreqFlanks.m3PFreqU[pIdx]);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_di_freq(unsigned pIdx, TM1MRNADiFreq<TRNAString>& pDiFreqs)
{
    mScaledFeatures[pIdx][9] = scale_feat_val(9, pDiFreqs.mSeedFreqUC[pIdx]);
    mScaledFeatures[pIdx][10] = scale_feat_val(10, pDiFreqs.mSeedFreqCA[pIdx]);
    mScaledFeatures[pIdx][11] = scale_feat_val(11, pDiFreqs.mSeedFreqCG[pIdx]);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_di_freq_flank(unsigned pIdx, TM1MRNADiFreqFlank<TRNAString>& pDiFreqFlanks)
{
    mScaledFeatures[pIdx][12] = scale_feat_val(12, pDiFreqFlanks.m3PFreqAA[pIdx]);
    mScaledFeatures[pIdx][13] = scale_feat_val(13, pDiFreqFlanks.m3PFreqAU[pIdx]);
    mScaledFeatures[pIdx][14] = scale_feat_val(14, pDiFreqFlanks.m3PFreqAG[pIdx]);
    mScaledFeatures[pIdx][15] = scale_feat_val(15, pDiFreqFlanks.m3PFreqAC[pIdx]);
    mScaledFeatures[pIdx][16] = scale_feat_val(16, pDiFreqFlanks.m3PFreqUA[pIdx]);
    mScaledFeatures[pIdx][17] = scale_feat_val(17, pDiFreqFlanks.m3PFreqUU[pIdx]);
    mScaledFeatures[pIdx][18] = scale_feat_val(18, pDiFreqFlanks.m3PFreqUG[pIdx]);
    mScaledFeatures[pIdx][19] = scale_feat_val(19, pDiFreqFlanks.m3PFreqUC[pIdx]);
    mScaledFeatures[pIdx][20] = scale_feat_val(20, pDiFreqFlanks.m3PFreqGU[pIdx]);
    mScaledFeatures[pIdx][21] = scale_feat_val(21, pDiFreqFlanks.m3PFreqCA[pIdx]);
    mScaledFeatures[pIdx][22] = scale_feat_val(22, pDiFreqFlanks.m3PFreqCU[pIdx]);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_single_match(unsigned pIdx, TM1MRNASingleMatch<TRNAString>& pSingleMatches)
{
    mScaledFeatures[pIdx][23] = scale_feat_val(23, pSingleMatches.mFreqUG[pIdx]);
    mScaledFeatures[pIdx][24] = scale_feat_val(24, pSingleMatches.mFreqGU[pIdx]);
    mScaledFeatures[pIdx][25] = scale_feat_val(25, pSingleMatches.mFreqCG[pIdx]);
}

template <class TRNAString>
void TM1ScaledFeatures<TRNAString>::scale_two_match(unsigned pIdx, TM1MRNATwoConsecMatch<TRNAString>& pTwoMatches)
{
    mScaledFeatures[pIdx][26] = scale_feat_val(26, pTwoMatches.mFreqUACG[pIdx]);
    mScaledFeatures[pIdx][27] = scale_feat_val(27, pTwoMatches.mFreqUAUG[pIdx]);
    mScaledFeatures[pIdx][28] = scale_feat_val(28, pTwoMatches.mFreqCGGC[pIdx]);
    mScaledFeatures[pIdx][29] = scale_feat_val(29, pTwoMatches.mFreqUGGC[pIdx]);
}

//
// TM1RawFeatures methods
//
template <class TRNAString>
void TM1MRNAFeatures<TRNAString>::clear_features()
{
    mSeedTypes.clear_features();
    mSiteCounts.clear_features();
    mAURich.clear_features();
    mSingleFreqs.clear_features();
    mSingleFreqFlanks.clear_features();
    mDiFreqs.clear_features();
    mDiFreqFlanks.clear_features();
    mSingleMatches.clear_features();
    mTwoMatches.clear_features();

    mScaledFeats.clear_features();
}

template <class TRNAString>
void TM1MRNAFeatures<TRNAString>::resize_features(unsigned pSize)
{
    mSeedTypes.resize_features(pSize);
    mSiteCounts.resize_features(pSize);
    mAURich.resize_features(pSize);
    mSingleFreqs.resize_features(pSize);
    mSingleFreqFlanks.resize_features(pSize);
    mDiFreqs.resize_features(pSize);
    mDiFreqFlanks.resize_features(pSize);
    mSingleMatches.resize_features(pSize);
    mTwoMatches.resize_features(pSize);

    mScaledFeats.resize_features(pSize);

}

template <class TRNAString>
int TM1MRNAFeatures<TRNAString>::add_features(
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures,
        TM1SortedSitePos<TRNAString> &pSortedSites)
{
    const StringSet<String<unsigned> >& sortedSites = pSortedSites.get_sorted_mrna_pos();
    const String<unsigned>& mRNAIDs = pSortedSites.get_mrna_ids();
    unsigned idx;

    resize_features((unsigned)length(mRNAIDs));

    for (unsigned i = 0; i < length(mRNAIDs); ++i)
    {
        idx = i;
        mSeedTypes.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mSiteCounts.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mAURich.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mSingleFreqs.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mSingleFreqFlanks.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mDiFreqs.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mDiFreqFlanks.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mSingleMatches.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
        mTwoMatches.add_features(idx, sortedSites[idx], pSeedSites, pRawFeatures);
    }

    mScaledFeats.scale_features(pSortedSites, mSeedTypes, mAURich, mSingleFreqs, mSingleFreqFlanks, mDiFreqs,
            mDiFreqFlanks, mSingleMatches, mTwoMatches);

//    print_features(pSortedSites);

    return 0;
}

template <class TRNAString>
void TM1MRNAFeatures<TRNAString>::print_features(TM1SortedSitePos<TRNAString> &pSortedSites)
{
    const String<unsigned>& mRNAIDs = pSortedSites.get_mrna_ids();

    for (unsigned i = 0; i < length(mRNAIDs); ++i)
    {
        mSiteCounts.print_feature(i);
        mSeedTypes.print_feature(i);
        mAURich.print_feature(i);
        mSingleFreqs.print_feature(i);
        mSingleFreqFlanks.print_feature(i);
        mDiFreqs.print_feature(i);
        mDiFreqFlanks.print_feature(i);
        mSingleMatches.print_feature(i);
        mTwoMatches.print_feature(i);
    }
    mScaledFeats.print_features();
    std::cout << std::endl;
}

//
// TM1MRNASeedType methods
//
template <class TRNAString>
void TM1MRNASeedType<TRNAString>::clear_features()
{
    clear(mNum6mer);
    clear(mNum8mer);
}

template <class TRNAString>
void TM1MRNASeedType<TRNAString>::resize_features(unsigned pSize)
{
    resize(mNum6mer, pSize);
    resize(mNum8mer, pSize);
    for (unsigned i = 0; i < pSize; ++i)
    {
        mNum6mer[i] = 0;
        mNum8mer[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNASeedType<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    CharString seedType;

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }
        seedType = pRawFeatures.get_seed_type(pSortedSites[i]);
        if (seedType == "6mer")
        {
            ++mNum6mer[pIdx];
        }
        else if (seedType == "8mer")
        {
            ++mNum8mer[pIdx];
        }
    }

    return 0;
}

template <class TRNAString>
void TM1MRNASeedType<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 1:" << mNum6mer[pIdx];
    std::cout << ", 10:" << mNum8mer[pIdx] << std::endl;
}

//
// TM1MRNASiteCount methods
//
template <class TRNAString>
void TM1MRNASiteCount<TRNAString>::clear_features()
{
    clear(mSiteCount);
}

template <class TRNAString>
void TM1MRNASiteCount<TRNAString>::resize_features(unsigned pSize)
{
    resize(mSiteCount, pSize);
    for (unsigned i = 0; i < pSize; ++i)
    {
        mSiteCount[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNASiteCount<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &)
{
    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }
        ++mSiteCount[pIdx];
    }

    return 0;
}

template <class TRNAString>
void TM1MRNASiteCount<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", site_count:" << mSiteCount[pIdx] << std::endl;
}

//
// TM1MRNAAURich methods
//
template <class TRNAString>
void TM1MRNAAURich<TRNAString>::clear_features()
{
    clear(mAU6mer);
    clear(mAU7merm8);
    clear(mAU8mer);
}

template <class TRNAString>
void TM1MRNAAURich<TRNAString>::resize_features(unsigned pSize)
{
    resize(mAU6mer, pSize);
    resize(mAU7merm8, pSize);
    resize(mAU8mer, pSize);

    for (unsigned i = 0; i < pSize; ++i)
    {
        mAU6mer[i] = 0;
        mAU7merm8[i] = 0;
        mAU8mer[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNAAURich<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    CharString seedType;
    const TM1FeatAURich<TRNAString>& auRich = pRawFeatures.get_au_rich();

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }

        if (auRich.mAURich[pSortedSites[i]] < 0.6)
        {
            continue;
        }

        seedType = pRawFeatures.get_seed_type(pSortedSites[i]);
        if (seedType == "6mer")
        {
            ++mAU6mer[pIdx];
        }
        else if (seedType == "7mer-m8")
        {
            ++mAU7merm8[pIdx];
        }
        else if (seedType == "8mer")
        {
            ++mAU8mer[pIdx];
        }
    }

    return 0;
}

template <class TRNAString>
void TM1MRNAAURich<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 2:" << mAU6mer[pIdx];
    std::cout << ", 5:" << mAU7merm8[pIdx];
    std::cout << ", 11:" << mAU8mer[pIdx] << std::endl;
}

//
// TM1MRNASingleFreq methods
//
template <class TRNAString>
void TM1MRNASingleFreq<TRNAString>::clear_features()
{
    clear(mSeedFreqU);
    clear(mSeedFreqC);
}

template <class TRNAString>
void TM1MRNASingleFreq<TRNAString>::resize_features(unsigned pSize)
{
    resize(mSeedFreqU, pSize);
    resize(mSeedFreqC, pSize);

    for (unsigned i = 0; i < pSize; ++i)
    {
        mSeedFreqU[i] = 0;
        mSeedFreqC[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNASingleFreq<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    const TM1FeatSingleFreq<TRNAString>& singleFreq = pRawFeatures.get_single_freq();
    CharString seedType;
    unsigned idx;
    float totalFreqU = 0;
    float totalFreqC = 0;
    float seedlen = 0;
    int site_count = 0;

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }
        ++site_count;
    }

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        idx = (unsigned)length(pSortedSites) - 1 - i;

        if (!pSeedSites.mEffectiveSites[pSortedSites[idx]])
        {
            continue;
        }

        totalFreqU += singleFreq.mSeedCountU[pSortedSites[idx]];
        totalFreqC += singleFreq.mSeedCountC[pSortedSites[idx]];

        seedType = pRawFeatures.get_seed_type(pSortedSites[idx]);
        if (seedType == "6mer")
        {
            seedlen = (float)(6 * site_count);
        }
        else if (seedType == "7mer-A1" || seedType == "7mer-m8")
        {
            seedlen = (float)(7 * site_count);
        }
        else if (seedType == "8mer")
        {
            seedlen =(float)(8 * site_count);
        }
        mSeedFreqU[pIdx] += totalFreqU / seedlen;
        mSeedFreqC[pIdx] += totalFreqC / seedlen;
    }

    return 0;
}

template <class TRNAString>
void TM1MRNASingleFreq<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 14:" << mSeedFreqU[pIdx];
    std::cout << ", 16:" << mSeedFreqC[pIdx] << std::endl;
}

//
// TM1MRNASingleFreqFlank methods
//
template <class TRNAString>
void TM1MRNASingleFreqFlank<TRNAString>::clear_features()
{
    clear(m3PFreqA);
    clear(m3PFreqU);
}

template <class TRNAString>
void TM1MRNASingleFreqFlank<TRNAString>::resize_features(unsigned pSize)
{
    resize(m3PFreqA, pSize);
    resize(m3PFreqU, pSize);

    for (unsigned i = 0; i < pSize; ++i)
    {
        m3PFreqA[i] = 0;
        m3PFreqU[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNASingleFreqFlank<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    const TM1FeatSingleFreqFlank<TRNAString>& singleFreqFlank = pRawFeatures.get_single_freq_flank();
    unsigned idx;
    float totalFreqA = 0;
    float totalFreqU = 0;
    float seqlen = 0;
    int site_count = 0;

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }
        ++site_count;
    }

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        idx = (unsigned)length(pSortedSites) - 1 - i;

        if (!pSeedSites.mEffectiveSites[pSortedSites[idx]])
        {
            continue;
        }

        totalFreqA += singleFreqFlank.m3PCountA[pSortedSites[idx]];
        totalFreqU += singleFreqFlank.m3PCountU[pSortedSites[idx]];

        seqlen = (float)(singleFreqFlank.m3PTotal[pSortedSites[idx]] * site_count);
        if (seqlen != 0)
        {
            m3PFreqA[pIdx] += totalFreqA / seqlen;
            m3PFreqU[pIdx] += totalFreqU / seqlen;
        }
    }

    return 0;
}

template <class TRNAString>
void TM1MRNASingleFreqFlank<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 17:" << m3PFreqA[pIdx];
    std::cout << ", 18:" << m3PFreqU[pIdx] << std::endl;
}

//
// TM1MRNADiFreq methods
//
template <class TRNAString>
void TM1MRNADiFreq<TRNAString>::clear_features()
{
    clear(mSeedFreqUC);
    clear(mSeedFreqCA);
    clear(mSeedFreqCG);
}

template <class TRNAString>
void TM1MRNADiFreq<TRNAString>::resize_features(unsigned pSize)
{
    resize(mSeedFreqUC, pSize);
    resize(mSeedFreqCA, pSize);
    resize(mSeedFreqCG, pSize);

    for (unsigned i = 0; i < pSize; ++i)
    {
        mSeedFreqUC[i] = 0;
        mSeedFreqCA[i] = 0;
        mSeedFreqCG[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNADiFreq<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    const TM1FeatDiFreq<TRNAString>& diFreq = pRawFeatures.get_di_freq();
    CharString seedType;
    unsigned idx;
    float totalFreqUC = 0;
    float totalFreqCA = 0;
    float totalFreqCG = 0;
    float seedlen = 0;
    int site_count = 0;

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }
        ++site_count;
    }


    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        idx = (unsigned)length(pSortedSites) - 1 - i;

        if (!pSeedSites.mEffectiveSites[pSortedSites[idx]])
        {
            continue;
        }

        totalFreqUC += diFreq.mSeedCountUC[pSortedSites[idx]];
        totalFreqCA += diFreq.mSeedCountCA[pSortedSites[idx]];
        totalFreqCG += diFreq.mSeedCountCG[pSortedSites[idx]];

        seedType = pRawFeatures.get_seed_type(pSortedSites[idx]);
        if (seedType == "6mer")
        {
            seedlen = (float)(6 * site_count);
        }
        else if (seedType == "7mer-A1" || seedType == "7mer-m8")
        {
            seedlen = (float)(7 * site_count);
        }
        else if (seedType == "8mer")
        {
            seedlen = (float)(8 * site_count);
        }
        mSeedFreqUC[pIdx] += totalFreqUC / seedlen;
        mSeedFreqCA[pIdx] += totalFreqCA / seedlen;
        mSeedFreqCG[pIdx] += totalFreqCG / seedlen;
    }

    return 0;
}

template <class TRNAString>
void TM1MRNADiFreq<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 28:" << mSeedFreqUC[pIdx];
    std::cout << ", 33:" << mSeedFreqCA[pIdx];
    std::cout << ", 35:" << mSeedFreqCG[pIdx] << std::endl;
}

//
// TM1MRNADiFreqFlank methods
//
template <class TRNAString>
void TM1MRNADiFreqFlank<TRNAString>::clear_features()
{
    clear(m3PFreqAA);
    clear(m3PFreqAU);
    clear(m3PFreqAG);
    clear(m3PFreqAC);
    clear(m3PFreqUA);
    clear(m3PFreqUU);
    clear(m3PFreqUG);
    clear(m3PFreqUC);
    clear(m3PFreqGU);
    clear(m3PFreqCA);
    clear(m3PFreqCU);
}

template <class TRNAString>
void TM1MRNADiFreqFlank<TRNAString>::resize_features(unsigned pSize)
{
    resize(m3PFreqAA, pSize);
    resize(m3PFreqAU, pSize);
    resize(m3PFreqAG, pSize);
    resize(m3PFreqAC, pSize);
    resize(m3PFreqUA, pSize);
    resize(m3PFreqUU, pSize);
    resize(m3PFreqUG, pSize);
    resize(m3PFreqUC, pSize);
    resize(m3PFreqGU, pSize);
    resize(m3PFreqCA, pSize);
    resize(m3PFreqCU, pSize);

    for (unsigned i = 0; i < pSize; ++i)
    {
        m3PFreqAA[i] = 0;
        m3PFreqAU[i] = 0;
        m3PFreqAG[i] = 0;
        m3PFreqAC[i] = 0;
        m3PFreqUA[i] = 0;
        m3PFreqUU[i] = 0;
        m3PFreqUG[i] = 0;
        m3PFreqUC[i] = 0;
        m3PFreqGU[i] = 0;
        m3PFreqCA[i] = 0;
        m3PFreqCU[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNADiFreqFlank<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    const TM1FeatDiFreqFlank<TRNAString>& diFreqFlank = pRawFeatures.get_di_freq_flank();
    unsigned idx;
    float totalFreqAA = 0;
    float totalFreqAU = 0;
    float totalFreqAG = 0;
    float totalFreqAC = 0;
    float totalFreqUA = 0;
    float totalFreqUU = 0;
    float totalFreqUG = 0;
    float totalFreqUC = 0;
    float totalFreqGU = 0;
    float totalFreqCA = 0;
    float totalFreqCU = 0;
    float seqlen = 0;
    int site_count = 0;

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }
        ++site_count;
    }

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        idx = (unsigned)length(pSortedSites) - 1 - i;

        if (!pSeedSites.mEffectiveSites[pSortedSites[idx]])
        {
            continue;
        }

        totalFreqAA += diFreqFlank.m3PCountAA[pSortedSites[idx]];
        totalFreqAU += diFreqFlank.m3PCountAU[pSortedSites[idx]];
        totalFreqAG += diFreqFlank.m3PCountAG[pSortedSites[idx]];
        totalFreqAC += diFreqFlank.m3PCountAC[pSortedSites[idx]];
        totalFreqUA += diFreqFlank.m3PCountUA[pSortedSites[idx]];
        totalFreqUU += diFreqFlank.m3PCountUU[pSortedSites[idx]];
        totalFreqUG += diFreqFlank.m3PCountUG[pSortedSites[idx]];
        totalFreqUC += diFreqFlank.m3PCountUC[pSortedSites[idx]];
        totalFreqGU += diFreqFlank.m3PCountGU[pSortedSites[idx]];
        totalFreqCA += diFreqFlank.m3PCountCA[pSortedSites[idx]];
        totalFreqCU += diFreqFlank.m3PCountCU[pSortedSites[idx]];

        seqlen = (float)(diFreqFlank.m3PTotal[pSortedSites[idx]] * site_count);
        if (seqlen != 0)
        {
            m3PFreqAA[pIdx] += totalFreqAA / seqlen;
            m3PFreqAU[pIdx] += totalFreqAU / seqlen;
            m3PFreqAG[pIdx] += totalFreqAG / seqlen;
            m3PFreqAC[pIdx] += totalFreqAC / seqlen;
            m3PFreqUA[pIdx] += totalFreqUA / seqlen;
            m3PFreqUU[pIdx] += totalFreqUU / seqlen;
            m3PFreqUG[pIdx] += totalFreqUG / seqlen;
            m3PFreqUC[pIdx] += totalFreqUC / seqlen;
            m3PFreqGU[pIdx] += totalFreqGU / seqlen;
            m3PFreqCA[pIdx] += totalFreqCA / seqlen;
            m3PFreqCU[pIdx] += totalFreqCU / seqlen;
        }

    }

    return 0;
}

template <class TRNAString>
void TM1MRNADiFreqFlank<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 37:" << m3PFreqAA[pIdx];
    std::cout << ", 38:" << m3PFreqAU[pIdx];
    std::cout << ", 39:" << m3PFreqAG[pIdx];
    std::cout << ", 40:" << m3PFreqAC[pIdx];
    std::cout << ", 41:" << m3PFreqUA[pIdx];
    std::cout << ", 42:" << m3PFreqUU[pIdx];
    std::cout << ", 43:" << m3PFreqUG[pIdx];
    std::cout << ", 44:" << m3PFreqUC[pIdx];
    std::cout << ", 46:" << m3PFreqGU[pIdx];
    std::cout << ", 49:" << m3PFreqCA[pIdx];
    std::cout << ", 50:" << m3PFreqCU[pIdx] << std::endl;
}

//
// TM1MRNASingleMatch methods
//
template <class TRNAString>
void TM1MRNASingleMatch<TRNAString>::clear_features()
{
    clear(mFreqUG);
    clear(mFreqGU);
    clear(mFreqCG);
}

template <class TRNAString>
void TM1MRNASingleMatch<TRNAString>::resize_features(unsigned pSize)
{
    resize(mFreqUG, pSize);
    resize(mFreqGU, pSize);
    resize(mFreqCG, pSize);

    for (unsigned i = 0; i < pSize; ++i)
    {
        mFreqUG[i] = 0;
        mFreqGU[i] = 0;
        mFreqCG[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNASingleMatch<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    const TM1FeatSingleMatch<TRNAString>& singleMatch = pRawFeatures.get_single_match();
    int totalUG = 0;
    int totalGU = 0;
    int totalCG = 0;
    int seqcount = 0;

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }

        totalUG += singleMatch.mMatchUG[pSortedSites[i]];
        totalGU += singleMatch.mMatchGU[pSortedSites[i]];
        totalCG += singleMatch.mMatchCG[pSortedSites[i]];

        ++seqcount;
    }

    mFreqUG[pIdx] = (float) totalUG / seqcount;
    mFreqGU[pIdx] = (float) totalGU / seqcount;
    mFreqCG[pIdx] = (float) totalCG / seqcount;

    return 0;
}

template <class TRNAString>
void TM1MRNASingleMatch<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 55:" << mFreqUG[pIdx];
    std::cout << ", 57:" << mFreqGU[pIdx];
    std::cout << ", 58:" << mFreqCG[pIdx] << std::endl;
}

//
// TM1MRNATwoConsecMatch methods
//
template <class TRNAString>
void TM1MRNATwoConsecMatch<TRNAString>::clear_features()
{
    clear(mFreqUACG);
    clear(mFreqUAUG);
    clear(mFreqCGGC);
    clear(mFreqUGGC);
}

template <class TRNAString>
void TM1MRNATwoConsecMatch<TRNAString>::resize_features(unsigned pSize)
{
    resize(mFreqUACG, pSize);
    resize(mFreqUAUG, pSize);
    resize(mFreqCGGC, pSize);
    resize(mFreqUGGC, pSize);

    for (unsigned i = 0; i < pSize; ++i)
    {
        mFreqUACG[i] = 0;
        mFreqUAUG[i] = 0;
        mFreqCGGC[i] = 0;
        mFreqUGGC[i] = 0;
    }
}

template <class TRNAString>
int TM1MRNATwoConsecMatch<TRNAString>::add_features(
        unsigned pIdx,
        const String<unsigned>& pSortedSites,
        TM1SeedSites<TRNAString> &pSeedSites,
        TM1RawFeatures<TRNAString> &pRawFeatures)
{
    const TM1FeatTwoConsecMatch<TRNAString>& twoConsecMatch = pRawFeatures.get_two_consec_match();
    int totalUACG = 0;
    int totalUAUG = 0;
    int totalCGGC = 0;
    int totalUGGC = 0;
    int seqcount = 0;

    for (unsigned i = 0; i < length(pSortedSites); ++i)
    {
        if (!pSeedSites.mEffectiveSites[pSortedSites[i]])
        {
            continue;
        }

        totalUACG += twoConsecMatch.mMatchUACG[pSortedSites[i]];
        totalUAUG += twoConsecMatch.mMatchUAUG[pSortedSites[i]];
        totalCGGC += twoConsecMatch.mMatchCGGC[pSortedSites[i]];
        totalUGGC += twoConsecMatch.mMatchUGGC[pSortedSites[i]];

        ++seqcount;
    }

    mFreqUACG[pIdx] = (float) totalUACG / seqcount;
    mFreqUAUG[pIdx] = (float) totalUAUG / seqcount;
    mFreqCGGC[pIdx] = (float) totalCGGC / seqcount;
    mFreqUGGC[pIdx] = (float) totalUGGC / seqcount;

    return 0;
}

template <class TRNAString>
void TM1MRNATwoConsecMatch<TRNAString>::print_feature(unsigned pIdx)
{
    std::cout << pIdx;
    std::cout << ", 68:" << mFreqUACG[pIdx];
    std::cout << ", 70:" << mFreqUAUG[pIdx];
    std::cout << ", 79:" << mFreqCGGC[pIdx];
    std::cout << ", 89:" << mFreqUGGC[pIdx] << std::endl;
}

// Explicit template instantiation
template class TM1MRNASeedType<TRNATYPE>;
template class TM1MRNASiteCount<TRNATYPE>;
template class TM1MRNAAURich<TRNATYPE>;
template class TM1MRNASingleFreq<TRNATYPE>;
template class TM1MRNASingleFreqFlank<TRNATYPE>;
template class TM1MRNADiFreq<TRNATYPE>;
template class TM1MRNADiFreqFlank<TRNATYPE>;
template class TM1MRNASingleMatch<TRNATYPE>;
template class TM1MRNATwoConsecMatch<TRNATYPE>;

template class TM1ScaledFeatures<TRNATYPE>;
template class TM1MRNAFeatures<TRNATYPE>;

} // namespace tm1p
