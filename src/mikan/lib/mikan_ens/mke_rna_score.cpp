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

    for (unsigned i = 0; i < length(scoreTypes); ++i) {
        seqan::CharString sType;
        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);
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

    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();
    const mikan::TMRNAPosSet &RNAPos = pRNAScores.get_mrna_pos();
    const mikan::TCharSet &scoreTypes = pRNAScores.get_score_types();

    for (unsigned i = 0; i < length(scoreTypes); ++i) {
        seqan::CharString sType;
        append(sType, pPrefix);
        append(sType, ":");
        append(sType, scoreTypes[i]);

        unsigned idxTool = mIdxMap[std::string(toCString(sType))];

        for (unsigned j = 0; j < length(pRNAScores.mEffectiveRNAs); j++) {
            if (!pRNAScores.mEffectiveRNAs[j]) {
                continue;
            }

            unsigned idxSite = pRNAWithSites.get_idx_from_pos(RNAPos[j]);
            float score = pRNAScores.get_score(i, j);
            mRNARawScoreList[idxTool][idxSite] = score;
            mRNANormScoreList[idxTool][idxSite] = normalize_score(score, pMKEOpts, sType);
            mMRNAPos[idxSite] = RNAPos[j];
            mEffectiveRNAs[idxSite] = true;
        }
    }

}

float MKERNAScores::normalize_score(
        float pScore,
        MKEOptions const &pMKEOpts,
        seqan::CharString &pScoreType) {

    float nScore = pScore;

    return nScore;
}

void MKERNAScores::combine_scores(MKEOptions const &pMKEOpts) {
    for (unsigned i = 0; i < length(mEffectiveRNAs); i++) {
        if (!mEffectiveRNAs[i]) {
            mToolScores[i] = "";
            mRNAScores[i] = 0;
            continue;
        }

        std::stringstream stream;
        float score = 0;
        for (unsigned j = 0; j < mScoreTypeN; ++j) {
            stream << mScoreTypes[j] << ":" ;
            float tscore = mRNARawScoreList[j][i];
            tscore = roundf(tscore * 100.0f) / 100.0f;
            stream << tscore << "," ;

            float weight = 1;
            score += weight * mRNANormScoreList[j][i];
        }

        mToolScores[i] = stream.str();
        mRNAScores[i] = score;

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
