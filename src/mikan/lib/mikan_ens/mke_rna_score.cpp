#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mke_rna_score.hpp"     // MKERNAScores

using namespace seqan;

namespace mkens {

//
// MKERNAScores methods
//
void MKERNAScores::clear_scores() {
    mikan::MKRNAScores::clear_scores();

    clear(mToolScores);

    for (unsigned i = 0; i < length(mRNARawScoreList); ++i) {
        clear(mRNARawScoreList[i]);
        clear(mRNANormScoreList[i]);
    }
}

void MKERNAScores::add_score_types(
        mkens::MKEOptions const &pMKEOpts,
        mikan::MKRNAScores &pRNAScores,
        seqan::CharString &pPrefix) {

    const mikan::TCharSet &scoreTypes = pRNAScores.get_score_types();
    const MKEConfig &conf = pMKEOpts.get_conf();

    for (unsigned i = 0; i < length(scoreTypes); ++i) {
        seqan::CharString sType, sTypeC;
        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);

        sTypeC = sType;
        replace(sType, 2, 3, ":rna:");
        std::string ckey = toCString(sTypeC);
        if (!conf.get_rna_flag(ckey)) {
            continue;
        }

        appendValue(mScoreTypes, sType);
        mIdxMap[std::string(toCString(mScoreTypes[mScoreTypeN]))] = mScoreTypeN;
        ++mScoreTypeN;
    }
}

void MKERNAScores::init_score_list(mkens::MKERMAWithSites &pRNAWithSites) {

    resize(mRNARawScoreList, mScoreTypeN);
    resize(mRNANormScoreList, mScoreTypeN);

    resize(mEffectiveRNAs, length(pRNAWithSites.mEffectiveRNAs), false);
    resize(mMRNAPos, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mRNAScores, length(pRNAWithSites.mEffectiveRNAs));
    resize(mSiteNum, length(pRNAWithSites.mEffectiveRNAs), 0);
    resize(mToolScores, length(pRNAWithSites.mEffectiveRNAs));
    for (unsigned i = 0; i < mScoreTypeN; ++i) {
        resize(mRNARawScoreList[i], length(pRNAWithSites.mEffectiveRNAs));
        resize(mRNANormScoreList[i], length(pRNAWithSites.mEffectiveRNAs));
        for (unsigned j = 0; j < length(pRNAWithSites.mEffectiveRNAs); ++j) {
            mRNARawScoreList[i][j] = 0;
            mRNANormScoreList[i][j] = 0;
        }
    }

}

void MKERNAScores::add_scores(
        MKEOptions const &pMKEOpts,
        mkens::MKERMAWithSites &pRNAWithSites,
        mikan::MKRNAScores &pRNAScores,
        seqan::CharString &pPrefix) {

    const mikan::TMRNAPosSet &RNAPos = pRNAScores.get_mrna_pos();
    const mikan::TCharSet &scoreTypes = pRNAScores.get_score_types();
    const MKEConfig &conf = pMKEOpts.get_conf();

    for (unsigned i = 0; i < length(scoreTypes); ++i) {
        seqan::CharString sType, sTypeC;
        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);

        sTypeC = sType;
        replace(sType, 2, 3, ":rna:");
        std::string ckey = toCString(sTypeC);
        if (!conf.get_rna_flag(ckey)) {
            continue;
        }
        float lBound = conf.get_rna_lower(ckey);
        float uBound = conf.get_rna_upper(ckey);
        bool isRev = conf.get_rna_reverse(ckey);

        unsigned idxTool = mIdxMap[std::string(toCString(sType))];

        for (unsigned j = 0; j < length(pRNAScores.mEffectiveRNAs); j++) {
            if (!pRNAScores.mEffectiveRNAs[j]) {
                continue;
            }

            unsigned idxSite = pRNAWithSites.get_idx_from_pos(RNAPos[j]);
            float score = pRNAScores.get_score(i, j);
            mRNARawScoreList[idxTool][idxSite] = score;
            mRNANormScoreList[idxTool][idxSite] = normalize_score(score, lBound, uBound, isRev);
            mMRNAPos[idxSite] = RNAPos[j];
            mEffectiveRNAs[idxSite] = true;
        }
    }

}

float MKERNAScores::normalize_score(
        float pScore,
        float pLower,
        float pUpper,
        bool pReverse) {

    float nScore;

    if (pReverse) {
        nScore = (pScore - pUpper) / (pLower - pUpper);

    } else {
        nScore = (pScore - pLower) / (pUpper - pLower);
    }

    if (nScore < 0) {
        nScore = 0;
    } else if (nScore > 1) {
        nScore = 1;
    }

    return nScore;
}

void MKERNAScores::combine_scores(MKEOptions const &pMKEOpts) {
    const MKEConfig &conf = pMKEOpts.get_conf();
    seqan::StringSet<float> weights;
    float total_weight = 0;

    resize(weights, mScoreTypeN);
    for (unsigned i = 0; i < mScoreTypeN; ++i) {
        seqan::CharString sType = mScoreTypes[i];
        replace(sType, 2, 3, ":rna:");
        std::string ckey = toCString(sType);
        weights[i] = conf.get_site_weight(ckey);
        total_weight += weights[i];
    }

    for (unsigned i = 0; i < length(mEffectiveRNAs); i++) {
        if (!mEffectiveRNAs[i]) {
            mToolScores[i] = "";
            mRNAScores[i] = 0;
            continue;
        }

        std::stringstream stream;
        float score = 0;
        for (unsigned j = 0; j < mScoreTypeN; ++j) {
            stream << mScoreTypes[j] << ":";
            float tscore = mRNARawScoreList[j][i];
            tscore = roundf(tscore * 100.0f) / 100.0f;
            stream << tscore << ",";

            score += weights[j] * mRNANormScoreList[j][i];
        }

        mToolScores[i] = stream.str();
        mRNAScores[i] = score / total_weight;

    }
}

void MKERNAScores::set_site_count(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &pSiteScores,
        mkens::MKERMAWithSites &pRNAWithSites) {

    const mikan::TMRNAPosSet &RNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(pSeedSites.mEffectiveSites); i++) {
        if (!pSeedSites.mEffectiveSites[i] || !pSiteScores.mEffectiveSites[i]) {
            continue;
        }
        unsigned idxSite = pRNAWithSites.get_idx_from_pos(RNAPos[i]);
        ++mSiteNum[idxSite];
    }

}

} // namespace mkens
