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
        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);

        sTypeC = sType;
        replace(sTypeC, 2, 3, mKeySep.c_str());
        std::string ckey = toCString(sTypeC);
        if (!conf.get_site_flag(ckey)) {
            continue;
        }

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
        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);
        sTypeC = sType;
        replace(sTypeC, 2, 3, mKeySep.c_str());

        std::string ckey = toCString(sTypeC);
        if (!conf.get_site_flag(ckey)) {
            continue;
        }
        float lBound = conf.get_site_lower(ckey);
        float uBound = conf.get_site_upper(ckey);
        bool isRev = conf.get_site_reverse(ckey);

        unsigned idxTool = mIdxMap[std::string(toCString(sType))];

        for (unsigned j = 0; j < length(pSeedSites.mEffectiveSites); j++) {
            if (!pSeedSites.mEffectiveSites[j] || !pSiteScores.mEffectiveSites[j]) {
                continue;
            }

            unsigned idxSite = pMKESeedSites.get_idx_from_pos(RNAPos[j], S1Pos[j]);
            float score = pSiteScores.get_score(i, j);
            mSiteRawScoreList[idxTool][idxSite] = score;
            mSiteNormScoreList[idxTool][idxSite] = normalize_score(score, lBound, uBound, isRev);
            mEffectiveSites[idxSite] = true;
        }
    }
}

float MKESiteScores::normalize_score(
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

void MKESiteScores::combine_scores(MKEOptions const &pMKEOpts) {
    const MKEConfig &conf = pMKEOpts.get_conf();
    seqan::StringSet<float> weights;
    float total_weight = 0;

    resize(weights, mScoreTypeN);
    for (unsigned i = 0; i < mScoreTypeN; ++i) {
        seqan::CharString sType = mScoreTypes[i];
        replace(sType, 2, 3, mKeySep.c_str());
        std::string ckey = toCString(sType);
        weights[i] = conf.get_site_weight(ckey);
        total_weight += weights[i];
    }

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
            tscore = roundf(tscore * 100.0f) / 100.0f;
            stream << tscore << ",";

            score += weights[j] * mSiteNormScoreList[j][i];
        }

        mToolScores[i] = stream.str();
        mSiteScores[i] = score / total_weight;

    }
}

void MKESiteScores::print_all_scores(MKEOptions const &pMKEOpts) {
    const MKEConfig &conf = pMKEOpts.get_conf();
    seqan::StringSet<float> weights;
    seqan::StringSet<float> lbounds;
    seqan::StringSet<float> ubounds;
    seqan::StringSet<bool> is_rev;
    float total_weight = 0;

    resize(weights, mScoreTypeN);
    resize(lbounds, mScoreTypeN);
    resize(ubounds, mScoreTypeN);
    resize(is_rev, mScoreTypeN);
    for (unsigned i = 0; i < mScoreTypeN; ++i) {
        seqan::CharString sType = mScoreTypes[i];
        replace(sType, 2, 3, mKeySep.c_str());
        std::string ckey = toCString(sType);
        weights[i] = conf.get_site_weight(ckey);
        lbounds[i] = conf.get_site_lower(ckey);
        ubounds[i] = conf.get_site_upper(ckey);
        is_rev[i] = conf.get_site_reverse(ckey);
        total_weight += weights[i];
    }

    std::cout << "### Site scores ###" << std::endl;
    std::cout << "total weight: " << total_weight << std::endl;

    for (unsigned i = 0; i < length(mEffectiveSites); i++) {
        if (!mEffectiveSites[i]) {
            continue;
        }
        std::cout << i << ". mk score: " << mSiteScores[i] << std::endl;

        std::cout << "\tnormalized score: " << std::endl;
        std::stringstream stream_n;
        bool first_n = true;
        for (unsigned j = 0; j < mScoreTypeN; ++j) {
            if (mSiteNormScoreList[j][i] == 0 && mSiteRawScoreList[j][i] == 0) {
                continue;
            }
            std::cout << "\t\t" << mScoreTypes[j] << ": ";
            std::cout << mSiteNormScoreList[j][i];
            std::cout << ", weight: " << weights[j] << std::endl;

            if (!first_n) {
                stream_n << " + ";
            }
            stream_n << "(" << mSiteNormScoreList[j][i] << " * " <<  weights[j] << ")";
            first_n = false;
        }
        std::cout << "\t\t(" << stream_n.str() << ") / " << total_weight << std::endl;

        std::cout << "\traw score: " << std::endl;
        for (unsigned j = 0; j < mScoreTypeN; ++j) {
            if (mSiteNormScoreList[j][i] == 0 && mSiteRawScoreList[j][i] == 0) {
                continue;
            }
            std::cout << "\t\t" << mScoreTypes[j] << ": ";
            std::cout << mSiteRawScoreList[j][i];
            std::cout << ", lower: " << lbounds[j];
            std::cout << ", uppper: " << ubounds[j] << std::endl;
            if (is_rev[j]) {
                std::cout << "\t\t(" << mSiteRawScoreList[j][i] << " - " << ubounds[j] << ") / ";
                std::cout << "(" << lbounds[j] << " - " << ubounds[j] << ")" << std::endl;
            } else {
                std::cout << "\t\t(" << mSiteRawScoreList[j][i] << " - " << lbounds[j] << ") / ";
                std::cout << "(" << ubounds[j] << " - " << lbounds[j] << ")" << std::endl;

            }
        }
    }

    std::cout << std::endl;

}

} // namespace mkens

