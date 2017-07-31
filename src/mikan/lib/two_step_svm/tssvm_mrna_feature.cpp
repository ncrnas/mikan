#include <math.h>                   // roundf
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_mrna_feature.hpp"   // TSSVMRNARawFeatures, TSSVMFeatUTRLen, TSSVMFeatSiteNum,

using namespace seqan;

namespace tssvm {

//
// TSSVMRawFeatures methods
//
void TSSVMRNARawFeatures::clear_features() {
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

void TSSVMRNARawFeatures::resize_feat(unsigned pLen) {
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

int TSSVMRNARawFeatures::add_features(
        TSSVMSeedSites &pSeedSites,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKRMAWithSites &pRNAWithSites,
        TSSVMSiteScores &pSiteScores) {

    mikan::TMRNAPosSet &mUniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();
    resize_feat(length(pRNAWithSites.mEffectiveRNAs));

    mikan::TSitePosSet sitePos;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }
            appendValue(sitePos, rnaSitePosMap[i][j]);
        }

        if (length(sitePos) == 0) {
            continue;
        }

        unsigned seqIdx = mUniqRNAPosSet[i];
        mEffectiveRNAs[i] = true;

        mUTRLen.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);
        mSiteNum.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);
        mTotDiscUTRLen.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);
        mSeedTypeNum.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);
        mDiscBin.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);
        mOptDist.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);
        mSiteNumFlg.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);
        mTotDisc.add_features(i, sitePos, pSeedSites, pMRNASeqs[seqIdx], pSiteScores);

        clear(sitePos);
    }

    return 0;
}

//
// TSSVMFeatUTRLen methods
//
void TSSVMFeatUTRLen::clear_features() {
    clear(mUTRLen);
}

int TSSVMFeatUTRLen::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &,
        TSSVMSeedSites &,
        mikan::TRNAStr const &pMRNASeq,
        TSSVMSiteScores &) {

    float mRNALen;

    resize(mUTRLen[pIdx], 1);
    mRNALen = (float) length(pMRNASeq);
    mRNALen = mRNALen / 22559;
    mUTRLen[pIdx][0] = roundf(mRNALen * 10000) / 10000;

    return 0;
}

//
// TSSVMFeatSiteNum methods
//
void TSSVMFeatSiteNum::clear_features() {
    clear(mSiteNum);
    clear(mSiteNumRaw);
}

int TSSVMFeatSiteNum::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &pSitePosSet,
        TSSVMSeedSites &,
        mikan::TRNAStr const &,
        TSSVMSiteScores &) {

    float siteNum;

    resize(mSiteNum[pIdx], 1);
    siteNum = (float) length(pSitePosSet);
    mSiteNumRaw[pIdx] = (int) siteNum;
    siteNum = siteNum / 38;
    mSiteNum[pIdx][0] = roundf(siteNum * 10000) / 10000;

    return 0;

}

//
// TSSVMFeatTotDiscUTRLen methods
//
void TSSVMFeatTotDiscUTRLen::clear_features() {
    clear(mTotDiscUTRLen);
}

int TSSVMFeatTotDiscUTRLen::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &pSitePosSet,
        TSSVMSeedSites &,
        mikan::TRNAStr const &pMRNASeq,
        TSSVMSiteScores &pSiteScores) {

    const seqan::String<float> &scores = pSiteScores.get_scores();
    float minVal = 0.0f;
    float shiftVal = -2.27441f;
    float totVal = 0.0f;
    float mRNALen;
    float discUTRLen;

    resize(mTotDiscUTRLen[pIdx], 1);

    mRNALen = (float) length(pMRNASeq);

    for (unsigned i = 0; i < length(pSitePosSet); ++i) {
        float score = scores[pSitePosSet[i]];
        score = score - shiftVal;
        if (score > minVal) {
            totVal += score;
        }
    }

    discUTRLen = totVal / mRNALen;
    discUTRLen = roundf(discUTRLen * 10000000) / 10000000;
    discUTRLen = roundf(discUTRLen * 1000000) / 1000000;

    mTotDiscUTRLen[pIdx] = discUTRLen;

    return 0;
}

//
// TSSVMFeatSeedTypeNum methods
//
void TSSVMFeatSeedTypeNum::clear_features() {
    clear(mSeedTypeNum);
}

int TSSVMFeatSeedTypeNum::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &pSitePosSet,
        TSSVMSeedSites &pSeedSites,
        mikan::TRNAStr const &,
        TSSVMSiteScores &) {

    float maxNum = 38;
    const StringSet<CharString> &seedTypes = pSeedSites.get_seed_types();

    resize(mSeedTypeNum[pIdx], 9, 0);

    for (unsigned i = 0; i < length(pSitePosSet); ++i) {

        if (seedTypes[pSitePosSet[i]] == "8mer") {
            mSeedTypeNum[pIdx][0] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "7mer-m8") {
            mSeedTypeNum[pIdx][1] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "7mer-A1") {
            mSeedTypeNum[pIdx][2] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "6mer") {
            mSeedTypeNum[pIdx][3] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "GUM") {
            mSeedTypeNum[pIdx][4] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "GUT") {
            mSeedTypeNum[pIdx][5] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "LP") {
            mSeedTypeNum[pIdx][6] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "BT") {
            mSeedTypeNum[pIdx][7] += 1;
        } else if (seedTypes[pSitePosSet[i]] == "BM") {
            mSeedTypeNum[pIdx][8] += 1;
        }
    }

    for (unsigned i = 0; i < length(mSeedTypeNum[pIdx]); ++i) {
        float seedTypeNum = mSeedTypeNum[pIdx][i];
        seedTypeNum /= maxNum;
        seedTypeNum = roundf(seedTypeNum * 10000) / 10000;
        mSeedTypeNum[pIdx][i] = seedTypeNum;
    }

    return 0;
}

//
// TSSVMFeatDiscBin methods
//
void TSSVMFeatDiscBin::clear_features() {
    clear(mDiscBin);
}

int TSSVMFeatDiscBin::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &pSitePosSet,
        TSSVMSeedSites &,
        mikan::TRNAStr const &,
        TSSVMSiteScores &pSiteScores) {

    const seqan::String<float> &scores = pSiteScores.get_scores();
    float maxNum = 38.0f;
    float shiftVal = -2.27441f;

    resize(mDiscBin[pIdx], 17, 0.0);

    for (unsigned i = 0; i < length(pSitePosSet); ++i) {
        float score = scores[pSitePosSet[i]];
        score = score - shiftVal;

        if (score > 3.4031112) {
            mDiscBin[pIdx][15] += 1;
        } else if (score > 3.3554676432) {
            mDiscBin[pIdx][14] += 1;
        } else if (score > 3.30748377528) {
            mDiscBin[pIdx][13] += 1;
        } else if (score > 3.25949990736) {
            mDiscBin[pIdx][12] += 1;
        } else if (score > 3.21151603944) {
            mDiscBin[pIdx][11] += 1;
        } else if (score > 3.14753754888) {
            mDiscBin[pIdx][10] += 1;
        } else if (score > 3.062459768888) {
            mDiscBin[pIdx][9] += 1;
        } else if (score > 2.9487958548) {
            mDiscBin[pIdx][8] += 1;
        } else if (score > 2.7973574064) {
            mDiscBin[pIdx][7] += 1;
        } else if (score > 2.59555291224) {
            mDiscBin[pIdx][6] += 1;
        } else if (score > 2.32636681632) {
            mDiscBin[pIdx][5] += 1;
        } else if (score > 1.96733858472) {
            mDiscBin[pIdx][4] += 1;
        } else if (score > 1.48886115) {
            mDiscBin[pIdx][3] += 1;
        } else if (score > 0.8507778) {
            mDiscBin[pIdx][2] += 1;
        } else if (score > 0) {
            mDiscBin[pIdx][1] += 1;
        } else {
            mDiscBin[pIdx][0] += 1;
        }
    }

    for (unsigned i = 0; i < length(mDiscBin[pIdx]); ++i) {
        float binCount = mDiscBin[pIdx][i];
        binCount /= maxNum;
        binCount = roundf(binCount * 1000000) / 1000000;
        mDiscBin[pIdx][i] = binCount;
    }

    return 0;
}

//
// TSSVMFeatOptDist methods
//
void TSSVMFeatOptDist::clear_features() {
    clear(mOptDist);
}

int TSSVMFeatOptDist::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &pSitePosSet,
        TSSVMSeedSites &pSeedSites,
        mikan::TRNAStr const &,
        TSSVMSiteScores &) {

    const String<unsigned> &mS8Pos = pSeedSites.get_site_pos_s8();
    unsigned prevPos, pos, posDiff;
    String<unsigned> dist;

    resize(mOptDist[pIdx], 2, 0.0);

    if (length(pSitePosSet) < 2) {
        return 0;
    }

    resize(dist, length(pSitePosSet));

    prevPos = (unsigned) mS8Pos[pSitePosSet[0]];
    for (unsigned i = 1; i < length(pSitePosSet); ++i) {
        pos = (unsigned) mS8Pos[pSitePosSet[i]];
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
            mOptDist[pIdx][0] += 1;
        }
        if (dist[i] >= 17 && dist[i] <= 25) {
            mOptDist[pIdx][1] += 1;
        }
    }

    mOptDist[pIdx][0] /= 100;
    mOptDist[pIdx][1] /= 100;

    return 0;
}

//
// TSSVMFeatSiteNumFlg methods
//
void TSSVMFeatSiteNumFlg::clear_features() {
    clear(mSiteNumFlg);
}

int TSSVMFeatSiteNumFlg::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &pSitePosSet,
        TSSVMSeedSites &,
        mikan::TRNAStr const &,
        TSSVMSiteScores &) {

    float siteNum;

    resize(mSiteNumFlg[pIdx], 3, 0.0);
    siteNum = (float) length(pSitePosSet);

    if (siteNum <= 1) {
        mSiteNumFlg[pIdx][0] = 1;
    } else if (siteNum >= 8) {
        mSiteNumFlg[pIdx][2] = 1;
    } else {
        mSiteNumFlg[pIdx][1] = 1;
    }

    return 0;
}

//
// TSSVMFeatTotDisc methods
//
void TSSVMFeatTotDisc::clear_features() {
    clear(mTotDisc);
}

int TSSVMFeatTotDisc::add_features(
        unsigned pIdx,
        mikan::TSitePosSet &pSitePosSet,
        TSSVMSeedSites &,
        mikan::TRNAStr const &,
        TSSVMSiteScores &pSiteScores) {

    const seqan::String<float> &scores = pSiteScores.get_scores();
    float divVal = 100.0f;
    float minVal = 0.0f;
    float shiftVal = -2.27441f;
    float totVal = 0.0f;

    resize(mTotDisc[pIdx], 1);

    for (unsigned i = 0; i < length(pSitePosSet); ++i) {
        float score = scores[pSitePosSet[i]];
        score = score - shiftVal;
        if (score > minVal) {
            totVal += score;
        }
    }

    totVal /= divVal;
    mTotDisc[pIdx] = totVal;

    return 0;
}

} // namespace tssvm
