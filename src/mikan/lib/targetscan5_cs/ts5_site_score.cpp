#include <math.h>                // roundf
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_feature.hpp"       // TS5RawFeatures
#include "ts5_site_score.hpp"    // TS5SiteScores, TS5ScoreSeedType, TS5FeatSitePos, TS5FeatAURich,

using namespace seqan;

namespace ts5cs {

//
// TS5SiteScores methods
//
void TS5SiteScores::resize_scores(int pSize) {
    resize(mSiteScores, pSize);

    mSeedTypes.resize_scores(pSize);
    mSitePos.resize_scores(pSize);
    mAURich.resize_scores(pSize);
    mThreePrimePair.resize_scores(pSize);
}

void TS5SiteScores::clear_scores() {
    mikan::MKSiteScores::clear_scores();

    mRawFeatures.clear_features();
    mSeedTypes.clear_scores();
    mSitePos.clear_scores();
    mAURich.clear_scores();
    mThreePrimePair.clear_scores();
}

int TS5SiteScores::calc_scores(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites &) {

    int retVal;
    mikan::TCharStr seedType;
    int sitePos;
    float auRich, threePrimePair;
    float seedTypeScore, sitePosScore, auRichScore, threePrimePairScore, totalScore;
    int lenScores;

    retVal = mRawFeatures.add_features(pMiRNASeq, pMRNASeqs, pSeedSites);
    if (retVal != 0) {
        return retVal;
    }

    const mikan::TCharSet &seedTypes = pSeedSites.get_seed_types();

    lenScores = (int) length(pSeedSites.mEffectiveSites);
    resize(mEffectiveSites, lenScores);
    resize_scores(lenScores);

    for (int i = 0; i < lenScores; ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            mSiteScores[i] = 0.0;
            mSeedTypes.set_score(i, 0.0);
            mSitePos.set_score(i, 0.0);
            mAURich.set_score(i, 0.0);
            mThreePrimePair.set_score(i, 0.0);
            mEffectiveSites[i] = false;
            continue;
        }

        mEffectiveSites[i] = true;

        seedType = seedTypes[i];
        sitePos = mRawFeatures.get_site_pos(i);
        auRich = mRawFeatures.get_au_rich(i);
        threePrimePair = mRawFeatures.get_three_prime_pair(i);

        seedTypeScore = mSeedTypes.calc_score(i, seedType);
        sitePosScore = mSitePos.calc_score(i, seedType, sitePos);
        auRichScore = mAURich.calc_score(i, seedType, auRich);
        threePrimePairScore = mThreePrimePair.calc_score(i, seedType, threePrimePair);

        totalScore = seedTypeScore + sitePosScore + auRichScore + threePrimePairScore;

        mSiteScores[i] = totalScore;

    }

    return 0;
}

//
// TS5ScoreSeedType methods
//

void TS5ScoreSeedType::resize_scores(int pSize) {
    resize(mSeedTypes, pSize);
}

void TS5ScoreSeedType::clear_scores() {
    clear(mSeedTypes);
}

float TS5ScoreSeedType::calc_score(int pIdx, mikan::TCharStr &pSeedType) {
    float score = 0.0;

    if (pSeedType == "8mer") {
        score = mContrib8mer;
    } else if (pSeedType == "7mer-A1") {
        score = mContrib7A1;
    } else if (pSeedType == "7mer-m8") {
        score = mContrib7M8;
    }

    mSeedTypes[pIdx] = score;

    return score;
}

//
// TS5FeatSitePos methods
//

void TS5ScoreSitePos::resize_scores(int pSize) {
    resize(mSitePos, pSize);
}

void TS5ScoreSitePos::clear_scores() {
    clear(mSitePos);
}

float TS5ScoreSitePos::calc_score(int pIdx, mikan::TCharStr &pSeedType, int &pSitePos) {
    float score = 0.0;

    if (pSeedType == "8mer") {
        score = mSlope8mer * (float) pSitePos + mICept8mer;
    } else if (pSeedType == "7mer-A1") {
        score = mSlope7A1 * (float) pSitePos + mICept7A1;
    } else if (pSeedType == "7mer-m8") {
        score = mSlope7M8 * (float) pSitePos + mICept7M8;
    }

    score = roundf(score * 1000.0f) / 1000.0f;
    mSitePos[pIdx] = score;

    return score;
}

//
// TS5FeatAURich methods
//

void TS5ScoreAURich::resize_scores(int pSize) {
    resize(mAURich, pSize);
}

void TS5ScoreAURich::clear_scores() {
    clear(mAURich);
}

float TS5ScoreAURich::calc_score(int pIdx, mikan::TCharStr &pSeedType, float &pAuRich) {
    float score = 0.0;

    if (pSeedType == "8mer") {
        score = mSlope8mer * pAuRich + mICept8mer;
    } else if (pSeedType == "7mer-A1") {
        score = mSlope7A1 * pAuRich + mICept7A1;
    } else if (pSeedType == "7mer-m8") {
        score = mSlope7M8 * pAuRich + mICept7M8;
    }

    score = roundf(score * 1000.0f) / 1000.0f;
    mAURich[pIdx] = score;

    return score;
}

//
// TS5FeatThreePrimePair methods
//

void TS5ScoreThreePrimePair::resize_scores(int pSize) {
    resize(mThreePrimePair, pSize);
}

void TS5ScoreThreePrimePair::clear_scores() {
    clear(mThreePrimePair);
}

float TS5ScoreThreePrimePair::calc_score(int pIdx, mikan::TCharStr &pSeedType, float &pThreePrimePair) {
    float score = 0.0;

    if (pSeedType == "8mer") {
        score = mSlope8mer * pThreePrimePair + mICept8mer;
    } else if (pSeedType == "7mer-A1") {
        score = mSlope7A1 * pThreePrimePair + mICept7A1;
    } else if (pSeedType == "7mer-m8") {
        score = mSlope7M8 * pThreePrimePair + mICept7M8;
    }

    mThreePrimePair[pIdx] = score;

    return score;
}

} // namespace ts5cs
