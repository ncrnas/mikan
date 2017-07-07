#include <math.h>                   // roundf
#include <mk_inst_template.hpp>     // TRNATYPE
#include <tssvm_mrna_feature.hpp>   // TSSVMRNARawFeatures, TSSVMFeatUTRLen, TSSVMFeatSiteNum,

using namespace seqan;
using namespace mikan;

namespace tssvm {

//
// TSSVMRawFeatures methods
//
template<class TRNAString>
void TSSVMRNARawFeatures<TRNAString>::clear_features() {
    clear(mEffectiveRNAs);
    mUTRLen.clear_features();
    mSiteNum.clear_features();
    mTotDiscUTRLen.clear_features();
    mSeedTypeNum.clear_features();
    mDiscBin.clear_features();
    mOptDist.clear_features();
    mSiteNumFlg.clear_features();
    mTotDisc.clear_features();

}

template<class TRNAString>
void TSSVMRNARawFeatures<TRNAString>::resize_feat(unsigned pLen) {
    resize(mEffectiveRNAs, pLen, false);
    mUTRLen.resize_feat(pLen);
    mSiteNum.resize_feat(pLen);
    mTotDiscUTRLen.resize_feat(pLen);
    mSeedTypeNum.resize_feat(pLen);
    mDiscBin.resize_feat(pLen);
    mOptDist.resize_feat(pLen);
    mSiteNumFlg.resize_feat(pLen);
    mTotDisc.resize_feat(pLen);
}

template<class TRNAString>
int TSSVMRNARawFeatures<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        StringSet<TRNAString> const &pMRNASeqs,
        TSSVMSeedSiteOverlap<TRNAString> &pOverlappedSites,
        TSSVMSiteInputVector<TRNAString> &pSiteInput) {
    TItSet itSet;
    std::set<unsigned> &rnaPosSet = pOverlappedSites.get_mrna_pos_set();
    StringSet<String<unsigned> > &sortedMRNAPos = pOverlappedSites.get_sorted_mrna_pos();
    String<unsigned> sitePosByMRNA;

    resize_feat(length(pMRNASeqs));

    for (itSet = rnaPosSet.begin(); itSet != rnaPosSet.end(); ++itSet) {
        clear(sitePosByMRNA);
        for (unsigned i = 0; i < length(sortedMRNAPos[*itSet]); ++i) {
            if (!pSeedSites.mEffectiveSites[sortedMRNAPos[*itSet][i]]) {
                continue;
            }
            appendValue(sitePosByMRNA, sortedMRNAPos[*itSet][i]);
        }

        if (length(sitePosByMRNA) == 0) {
            continue;
        }
        mEffectiveRNAs[*itSet] = true;

        mUTRLen.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
        mSiteNum.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
        mTotDiscUTRLen.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
        mSeedTypeNum.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
        mDiscBin.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
        mOptDist.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
        mSiteNumFlg.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
        mTotDisc.add_features(*itSet, sitePosByMRNA, pSeedSites, pMRNASeqs, pSiteInput);
    }

    return 0;
}

//
// TSSVMFeatUTRLen methods
//
template<class TRNAString>
void TSSVMFeatUTRLen<TRNAString>::clear_features() {
    clear(mUTRLen);
}

template<class TRNAString>
int TSSVMFeatUTRLen<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &,
        TSSVMSeedSites<TRNAString> &,
        StringSet<TRNAString> const &pMRNASeqs,
        TSSVMSiteInputVector<TRNAString> &) {
    float mRNALen;

    resize(mUTRLen[pMRNAPosIdx], 1);
    mRNALen = (float) length(pMRNASeqs[pMRNAPosIdx]);
    mRNALen = mRNALen / 22559;
    mUTRLen[pMRNAPosIdx][0] = roundf(mRNALen * 10000) / 10000;

    return 0;
}

//
// TSSVMFeatSiteNum methods
//
template<class TRNAString>
void TSSVMFeatSiteNum<TRNAString>::clear_features() {
    clear(mSiteNum);
    clear(mSiteNumRaw);
}

template<class TRNAString>
int TSSVMFeatSiteNum<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &pSitePosByMRNA,
        TSSVMSeedSites<TRNAString> &,
        StringSet<TRNAString> const &,
        TSSVMSiteInputVector<TRNAString> &) {
    float siteNum;

    resize(mSiteNum[pMRNAPosIdx], 1);
    siteNum = (float) length(pSitePosByMRNA);
    mSiteNumRaw[pMRNAPosIdx] = (int) siteNum;
    siteNum = siteNum / 38;
    mSiteNum[pMRNAPosIdx][0] = roundf(siteNum * 10000) / 10000;

    return 0;

}

//
// TSSVMFeatTotDiscUTRLen methods
//
template<class TRNAString>
void TSSVMFeatTotDiscUTRLen<TRNAString>::clear_features() {
    clear(mTotDiscUTRLen);
}

template<class TRNAString>
int TSSVMFeatTotDiscUTRLen<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &pSitePosByMRNA,
        TSSVMSeedSites<TRNAString> &,
        StringSet<TRNAString> const &pMRNASeqs,
        TSSVMSiteInputVector<TRNAString> &pSiteInput) {
    const seqan::String<float> &scors = pSiteInput.get_scores();
    float minVal = 0.0f;
    float shiftVal = -2.27441f;
    float totVal = 0.0f;
    float mRNALen;
    float discUTRLen;

    resize(mTotDiscUTRLen[pMRNAPosIdx], 1);

    mRNALen = (float) length(pMRNASeqs[pMRNAPosIdx]);

    for (unsigned i = 0; i < length(pSitePosByMRNA); ++i) {
        float score = scors[pSitePosByMRNA[i]];
        score = score - shiftVal;
        if (score > minVal) {
            totVal += score;
        }
    }

    discUTRLen = totVal / mRNALen;
    discUTRLen = roundf(discUTRLen * 10000000) / 10000000;
    discUTRLen = roundf(discUTRLen * 1000000) / 1000000;

    mTotDiscUTRLen[pMRNAPosIdx] = discUTRLen;

    return 0;
}

//
// TSSVMFeatSeedTypeNum methods
//
template<class TRNAString>
void TSSVMFeatSeedTypeNum<TRNAString>::clear_features() {
    clear(mSeedTypeNum);
}

template<class TRNAString>
int TSSVMFeatSeedTypeNum<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &pSitePosByMRNA,
        TSSVMSeedSites<TRNAString> &pSeedSites,
        StringSet<TRNAString> const &,
        TSSVMSiteInputVector<TRNAString> &) {
    float maxNum = 38;
    const StringSet<CharString> &seedTypes = pSeedSites.get_seed_types();

    resize(mSeedTypeNum[pMRNAPosIdx], 9, 0);

    for (unsigned i = 0; i < length(pSitePosByMRNA); ++i) {

        if (seedTypes[pSitePosByMRNA[i]] == "8mer") {
            mSeedTypeNum[pMRNAPosIdx][0] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "7mer-m8") {
            mSeedTypeNum[pMRNAPosIdx][1] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "7mer-A1") {
            mSeedTypeNum[pMRNAPosIdx][2] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "6mer") {
            mSeedTypeNum[pMRNAPosIdx][3] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "GUM") {
            mSeedTypeNum[pMRNAPosIdx][4] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "GUT") {
            mSeedTypeNum[pMRNAPosIdx][5] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "LP") {
            mSeedTypeNum[pMRNAPosIdx][6] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "BT") {
            mSeedTypeNum[pMRNAPosIdx][7] += 1;
        } else if (seedTypes[pSitePosByMRNA[i]] == "BM") {
            mSeedTypeNum[pMRNAPosIdx][8] += 1;
        }
    }

    for (unsigned i = 0; i < length(mSeedTypeNum[pMRNAPosIdx]); ++i) {
        float seedTypeNum = mSeedTypeNum[pMRNAPosIdx][i];
        seedTypeNum /= maxNum;
        seedTypeNum = roundf(seedTypeNum * 10000) / 10000;
        mSeedTypeNum[pMRNAPosIdx][i] = seedTypeNum;
    }

    return 0;
}

//
// TSSVMFeatDiscBin methods
//
template<class TRNAString>
void TSSVMFeatDiscBin<TRNAString>::clear_features() {
    clear(mDiscBin);
}

template<class TRNAString>
int TSSVMFeatDiscBin<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &pSitePosByMRNA,
        TSSVMSeedSites<TRNAString> &,
        StringSet<TRNAString> const &,
        TSSVMSiteInputVector<TRNAString> &pSiteInput) {
    const seqan::String<float> &scors = pSiteInput.get_scores();
    float maxNum = 38.0f;
    float shiftVal = -2.27441f;

    resize(mDiscBin[pMRNAPosIdx], 17, 0.0);

    for (unsigned i = 0; i < length(pSitePosByMRNA); ++i) {
        float score = scors[pSitePosByMRNA[i]];
        score = score - shiftVal;

        if (score > 3.4031112) {
            mDiscBin[pMRNAPosIdx][15] += 1;
        } else if (score > 3.3554676432) {
            mDiscBin[pMRNAPosIdx][14] += 1;
        } else if (score > 3.30748377528) {
            mDiscBin[pMRNAPosIdx][13] += 1;
        } else if (score > 3.25949990736) {
            mDiscBin[pMRNAPosIdx][12] += 1;
        } else if (score > 3.21151603944) {
            mDiscBin[pMRNAPosIdx][11] += 1;
        } else if (score > 3.14753754888) {
            mDiscBin[pMRNAPosIdx][10] += 1;
        } else if (score > 3.062459768888) {
            mDiscBin[pMRNAPosIdx][9] += 1;
        } else if (score > 2.9487958548) {
            mDiscBin[pMRNAPosIdx][8] += 1;
        } else if (score > 2.7973574064) {
            mDiscBin[pMRNAPosIdx][7] += 1;
        } else if (score > 2.59555291224) {
            mDiscBin[pMRNAPosIdx][6] += 1;
        } else if (score > 2.32636681632) {
            mDiscBin[pMRNAPosIdx][5] += 1;
        } else if (score > 1.96733858472) {
            mDiscBin[pMRNAPosIdx][4] += 1;
        } else if (score > 1.48886115) {
            mDiscBin[pMRNAPosIdx][3] += 1;
        } else if (score > 0.8507778) {
            mDiscBin[pMRNAPosIdx][2] += 1;
        } else if (score > 0) {
            mDiscBin[pMRNAPosIdx][1] += 1;
        } else {
            mDiscBin[pMRNAPosIdx][0] += 1;
        }
    }

    for (unsigned i = 0; i < length(mDiscBin[pMRNAPosIdx]); ++i) {
        float binCount = mDiscBin[pMRNAPosIdx][i];
        binCount /= maxNum;
        binCount = roundf(binCount * 1000000) / 1000000;
        mDiscBin[pMRNAPosIdx][i] = binCount;
    }

    return 0;
}

//
// TSSVMFeatOptDist methods
//
template<class TRNAString>
void TSSVMFeatOptDist<TRNAString>::clear_features() {
    clear(mOptDist);
}

template<class TRNAString>
int TSSVMFeatOptDist<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &pSitePosByMRNA,
        TSSVMSeedSites<TRNAString> &pSeedSites,
        StringSet<TRNAString> const &,
        TSSVMSiteInputVector<TRNAString> &) {
    const String<unsigned> &mS8Pos = pSeedSites.get_site_pos_s8();
    unsigned prevPos, pos, posDiff;
    String<unsigned> dist;

    resize(mOptDist[pMRNAPosIdx], 2, 0.0);

    if (length(pSitePosByMRNA) < 2) {
        return 0;
    }

    resize(dist, length(pSitePosByMRNA));

    prevPos = (unsigned) mS8Pos[pSitePosByMRNA[0]];
    for (unsigned i = 1; i < length(pSitePosByMRNA); ++i) {
        pos = (unsigned) mS8Pos[pSitePosByMRNA[i]];
        posDiff = pos - prevPos;
        dist[i] = posDiff;

        if (i == 1) {
            dist[0] = posDiff;
        } else if (posDiff < dist[i - 1]) {
            dist[i - 1] = posDiff;
        }

        prevPos = pos;
    }

    for (unsigned i = 0; i < length(dist); ++i) {

        if (dist[i] >= 14 && dist[i] <= 46) {
            mOptDist[pMRNAPosIdx][0] += 1;
        }
        if (dist[i] >= 17 && dist[i] <= 25) {
            mOptDist[pMRNAPosIdx][1] += 1;
        }
    }

    mOptDist[pMRNAPosIdx][0] /= 100;
    mOptDist[pMRNAPosIdx][1] /= 100;

    return 0;
}

//
// TSSVMFeatSiteNumFlg methods
//
template<class TRNAString>
void TSSVMFeatSiteNumFlg<TRNAString>::clear_features() {
    clear(mSiteNumFlg);
}

template<class TRNAString>
int TSSVMFeatSiteNumFlg<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &pSitePosByMRNA,
        TSSVMSeedSites<TRNAString> &,
        StringSet<TRNAString> const &,
        TSSVMSiteInputVector<TRNAString> &) {
    float siteNum;

    resize(mSiteNumFlg[pMRNAPosIdx], 3, 0.0);
    siteNum = (float) length(pSitePosByMRNA);

    if (siteNum <= 1) {
        mSiteNumFlg[pMRNAPosIdx][0] = 1;
    } else if (siteNum >= 8) {
        mSiteNumFlg[pMRNAPosIdx][2] = 1;
    } else {
        mSiteNumFlg[pMRNAPosIdx][1] = 1;
    }

    return 0;
}

//
// TSSVMFeatTotDisc methods
//
template<class TRNAString>
void TSSVMFeatTotDisc<TRNAString>::clear_features() {
    clear(mTotDisc);
}

template<class TRNAString>
int TSSVMFeatTotDisc<TRNAString>::add_features(
        unsigned pMRNAPosIdx,
        String<unsigned> &pSitePosByMRNA,
        TSSVMSeedSites<TRNAString> &,
        StringSet<TRNAString> const &,
        TSSVMSiteInputVector<TRNAString> &pSiteInput) {
    const seqan::String<float> &scors = pSiteInput.get_scores();
    float divVal = 100.0f;
    float minVal = 0.0f;
    float shiftVal = -2.27441f;
    float totVal = 0.0f;

    resize(mTotDisc[pMRNAPosIdx], 1);

    for (unsigned i = 0; i < length(pSitePosByMRNA); ++i) {
        float score = scors[pSitePosByMRNA[i]];
        score = score - shiftVal;
        if (score > minVal) {
            totVal += score;
        }
    }

    totVal /= divVal;
    mTotDisc[pMRNAPosIdx] = totVal;

    return 0;
}

// Explicit template instantiation
template
class TSSVMFeatUTRLen<TRNATYPE>;

template
class TSSVMFeatSiteNum<TRNATYPE>;

template
class TSSVMFeatTotDiscUTRLen<TRNATYPE>;

template
class TSSVMFeatSeedTypeNum<TRNATYPE>;

template
class TSSVMFeatDiscBin<TRNATYPE>;

template
class TSSVMFeatOptDist<TRNATYPE>;

template
class TSSVMFeatSiteNumFlg<TRNATYPE>;

template
class TSSVMFeatTotDisc<TRNATYPE>;

template
class TSSVMRNARawFeatures<TRNATYPE>;

} // namespace tssvm
