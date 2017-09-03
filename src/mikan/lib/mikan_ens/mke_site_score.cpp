#include <math.h>                // roundf
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mke_site_score.hpp"    // MKESiteScores

using namespace seqan;

namespace mkens {

//
// MKESiteScores methods
//
void MKESiteScores::clear_scores() {
    mikan::MKSiteScores::clear_scores();

    clear(mToolScores);

    for (unsigned i = 0; i < length(mSiteRawScoreList); ++i) {
        clear(mSiteRawScoreList[i]);
        clear(mSiteNormScoreList[i]);
    }
}

void MKESiteScores::add_score_types(
        MKEOptions const &pMKEOpts,
        mikan::MKSiteScores &pSiteScores,
        seqan::CharString &pPrefix) {

    const mikan::TCharSet &scoreTypes = pSiteScores.get_score_types();
    const MKEConfig &conf = pMKEOpts.get_conf();

    for (unsigned i = 0; i < length(scoreTypes); ++i) {
        seqan::CharString sType, sTypeC;
        append(sTypeC, pPrefix);
        append(sTypeC, ":site:");
        append(sTypeC, scoreTypes[i]);

        std::string ckey = toCString(sTypeC);
        if (!conf.get_site_flag(ckey)) {
            continue;
        }

        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);
        appendValue(mScoreTypes, sType);

        mIdxMap[std::string(toCString(mScoreTypes[mScoreTypeN]))] = mScoreTypeN;
        ++mScoreTypeN;
    }
}

void MKESiteScores::init_score_list(MKESeedSites &pMKESeedSites) {
    resize(mSiteRawScoreList, mScoreTypeN);
    resize(mSiteNormScoreList, mScoreTypeN);

    resize(mEffectiveSites, length(pMKESeedSites.mEffectiveSites), false);
    resize(mSiteScores, length(pMKESeedSites.mEffectiveSites));
    resize(mToolScores, length(pMKESeedSites.mEffectiveSites));
    for (unsigned i = 0; i < mScoreTypeN; ++i) {
        resize(mSiteRawScoreList[i], length(pMKESeedSites.mEffectiveSites));
        resize(mSiteNormScoreList[i], length(pMKESeedSites.mEffectiveSites));
        for (unsigned j = 0; j < length(pMKESeedSites.mEffectiveSites); ++j) {
            mSiteRawScoreList[i][j] = 0;
            mSiteNormScoreList[i][j] = 0;
        }
    }

}

void MKESiteScores::add_scores(
        MKEOptions const &pMKEOpts,
        mikan::MKSeedSites &pSeedSites,
        MKESeedSites &pMKESeedSites,
        mikan::MKSiteScores &pSiteScores,
        seqan::CharString &pPrefix) {

    const mikan::TMRNAPosSet &RNAPos = pSeedSites.get_mrna_pos();
    const mikan::TMRNAPosSet &S1Pos = pSeedSites.get_site_pos_s1();
    const mikan::TCharSet &scoreTypes = pSiteScores.get_score_types();
    const MKEConfig &conf = pMKEOpts.get_conf();

    for (unsigned i = 0; i < length(scoreTypes); ++i) {
        seqan::CharString sType, sTypeC;
        append(sTypeC, pPrefix);
        append(sTypeC, ":site:");
        append(sTypeC, scoreTypes[i]);

        std::string ckey = toCString(sTypeC);
        if (!conf.get_site_flag(ckey)) {
            continue;
        }
        float lBound = conf.get_site_lower(ckey);
        float uBound = conf.get_site_upper(ckey);

        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);
        unsigned idxTool = mIdxMap[std::string(toCString(sType))];

        for (unsigned j = 0; j < length(pSeedSites.mEffectiveSites); j++) {
            if (!pSeedSites.mEffectiveSites[j] || !pSiteScores.mEffectiveSites[j]) {
                continue;
            }

            unsigned idxSite = pMKESeedSites.get_idx_from_pos(RNAPos[j], S1Pos[j]);
            float score = pSiteScores.get_score(i, j);
            mSiteRawScoreList[idxTool][idxSite] = score;
            mSiteNormScoreList[idxTool][idxSite] = normalize_score(score, pMKEOpts, sTypeC);
            mEffectiveSites[idxSite] = true;
        }
    }

}

float MKESiteScores::normalize_score(
        float pScore,
        MKEOptions const &pMKEOpts,
        seqan::CharString &pScoreType) {

    float nScore = pScore;
    const MKEConfig &conf = pMKEOpts.get_conf();

    return nScore;
}

void MKESiteScores::combine_scores(MKEOptions const &pMKEOpts) {
    for (unsigned i = 0; i < length(mEffectiveSites); i++) {
        if (!mEffectiveSites[i]) {
            mToolScores[i] = "";
            mSiteScores[i] = 0;
            continue;
        }

        std::stringstream stream;
        float score = 0;
        for (unsigned j = 0; j < mScoreTypeN; ++j) {
            stream << mScoreTypes[j] << ":";
            float tscore = mSiteRawScoreList[j][i];
//            std::cout << i << ", " << j <<  ", " << tscore << std::endl;
            tscore = roundf(tscore * 100.0f) / 100.0f;
            stream << tscore << ",";

            float weight = 1;
            score += weight * mSiteNormScoreList[j][i];

        }

        mToolScores[i] = stream.str();
        mSiteScores[i] = score;

    }
}


} // namespace mkens

