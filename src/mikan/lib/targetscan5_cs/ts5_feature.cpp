#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_feature.hpp"       // TS5RawFeatures, TS5FeatSeedType, TS5FeatSitePos, TS5FeatAURich,

using namespace seqan;

namespace ts5cs {

//
// TS5RawFeatures methods
//
void TS5RawFeatures::clear_features() {
    mSeedTypes.clear_features();
    mSitePos.clear_features();
    mAURich.clear_features();
    mThreePrimePair.clear_features();
}

int TS5RawFeatures::add_features(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        TS5SeedSites const &pSeedSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();

    resize(mEffectiveSites, length(mRNAPos));
    mSeedTypes.add_features(pMiRNASeq, pMRNASeqs, mEffectiveSites, mRNAPos, sitePos);
    mSitePos.add_features(pMiRNASeq, pMRNASeqs, mEffectiveSites, mRNAPos, sitePos, mSeedTypes);
    mAURich.add_features(pMiRNASeq, pMRNASeqs, mEffectiveSites, mRNAPos, sitePos, mSeedTypes);
    mThreePrimePair.add_features(pMiRNASeq, pMRNASeqs, mEffectiveSites, mRNAPos, sitePos, mSeedTypes);

    return 0;
}

//
// TS5FeatSeedType methods
//
void TS5FeatSeedType::clear_features() {
    clear(mSeedTypes);
}

int TS5FeatSeedType::add_features(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        mikan::TSitePosSet const &pMRNAPos,
        mikan::TSitePosSet const &pSitePos) {
    bool IsA1, MatchM8;
    int startPos;
    unsigned endPos;
    mikan::TRNAStr revMiRNASeq, miRNAM8, mRNAM8, mRNAA1;

    revMiRNASeq = pMiRNASeq;
    complement(revMiRNASeq);
    miRNAM8 = revMiRNASeq[7];

    resize(mSeedTypes, length(pMRNAPos));
    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        // check m8 match
        startPos = pSitePos[i];
        MatchM8 = false;
        if (startPos > 0) {
            mRNAM8 = pMRNASeqs[pMRNAPos[i]][startPos - 1];
            if (miRNAM8 == mRNAM8) {
                MatchM8 = true;
            }
        }

        // check A1 site
        IsA1 = false;
        endPos = (unsigned) startPos + 6;
        if (endPos < length(pMRNASeqs[pMRNAPos[i]])) {
            mRNAA1 = pMRNASeqs[pMRNAPos[i]][endPos];

            if (mRNAA1 == 'A') {
                IsA1 = true;
            }
        }

        pEffectiveSites[i] = true;
        if (IsA1 && MatchM8) {
            mSeedTypes[i] = "8mer";
        } else if (IsA1) {
            mSeedTypes[i] = "7mer-A1";
        } else if (MatchM8) {
            mSeedTypes[i] = "7mer-m8";
        } else {
            mSeedTypes[i] = "6mer";
            pEffectiveSites[i] = false;
        }
    }

    return 0;
}

//
// TS5FeatSitePos methods
//
void TS5FeatSitePos::clear_features() {
    clear(mSitePos);
}

int TS5FeatSitePos::add_features(
        mikan::TRNAStr const &,
        mikan::TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        mikan::TSitePosSet const &pMRNAPos,
        mikan::TSitePosSet const &pSitePos,
        TS5FeatSeedType &pSeedTypes) {
    int seqLen, lenUp, lenDown, scoreLen;

    resize(mSitePos, length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            mSitePos[i] = -1;
            continue;
        }

        // Alias
        CharString const &seedType = pSeedTypes.get_seed_type(i);

        // Get length to UTR start
        lenUp = pSitePos[i];
        if (seedType == "8mer" || seedType == "7mer-m8") {
            lenUp -= 1;
        }
        if (lenUp < (MIN_DIST_TO_CDS - 1)) {
            pEffectiveSites[i] = false;
            mSitePos[i] = -1;
            continue;
        }

        // Get length to UTR end
        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        lenDown = seqLen - (lenUp + 7);

        // Get score
        scoreLen = std::min(lenUp, lenDown);
        scoreLen = std::min(scoreLen, mMaxLen);

        mSitePos[i] = scoreLen;

    }

    return 0;
}

//
// TS5FeatAURich methods
//
void TS5FeatAURich::clear_features() {
    clear(mAURich);
}

int TS5FeatAURich::add_features(
        mikan::TRNAStr const &,
        mikan::TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        mikan::TSitePosSet const &pMRNAPos,
        mikan::TSitePosSet const &pSitePos,
        TS5FeatSeedType &pSeedTypes) {
    int seqLen, startU, endU, startD, endD;
    float totalScore, upTotalScore, upMaxScore, downTotalScore, downMaxScore;
    CharString chrUp = "up";
    CharString chrDown = "down";

    resize(mAURich, length(pMRNAPos));
    startU = 0;
    startD = 0;
    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            mAURich[i] = -1.0;
            continue;
        }

        // Alias
        CharString const &seedType = pSeedTypes.get_seed_type(i);

        // Get start and end positions for upstream and downstream
        getUpDownStreamPos(seedType, pSitePos[i], startU, startD);
        endU = pSitePos[i] - 1;
        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        endD = std::min(startD + 30, seqLen);

        calcPosScores(seedType, chrUp, pMRNASeqs[pMRNAPos[i]], startU, endU, upTotalScore, upMaxScore);
        calcPosScores(seedType, chrDown, pMRNASeqs[pMRNAPos[i]], startD, endD, downTotalScore, downMaxScore);

        totalScore = (upTotalScore + downTotalScore) / (upMaxScore + downMaxScore);

        mAURich[i] = totalScore;
    }

    return 0;
}

void TS5FeatAURich::getUpDownStreamPos(
        CharString pSeedType,
        int pStartPos,
        int &pStartU,
        int &pStartD) {
    pStartU = std::max(0, pStartPos - 31);

    if (pSeedType == "8mer") {
        pStartD = pStartPos + 6;
    } else if (pSeedType == "7mer-A1") {
        pStartD = pStartPos + 7;
    } else if (pSeedType == "7mer-m8") {
        pStartD = pStartPos + 6;
    }
}

void TS5FeatAURich::calcPosScores(
        const CharString &pSeedType,
        CharString &pUpOrDown,
        const mikan::TRNAStr &pMRNASeq,
        int pStart,
        int pEnd,
        float &pTotalScore,
        float &pMaxScore) {
    float score = 0.0;
    int seqPos;

    pTotalScore = 0.0;
    pMaxScore = 0.0;

    for (int i = 0; i < (pEnd - pStart); ++i) {
        if (pSeedType == "8mer") {
            if (pUpOrDown == "up") {
                score = 1.0f / (i + 1);
            } else {
                score = 1.0f / (i + 2);
            }
        } else if (pSeedType == "7mer-A1") {
            score = 1.0f / (i + 2);
        } else if (pSeedType == "7mer-m8") {
            if (pUpOrDown == "up") {
                score = 1.0f / (i + 1);
            } else {
                if (i == 0) {
                    score = 0.5;
                } else {
                    score = 1.0f / (i + 1);
                }
            }
        }

        if (pUpOrDown == "up") {
            seqPos = pEnd - i - 1;
        } else {
            seqPos = i + pStart;
        }

        if (pMRNASeq[seqPos] == 'A' || pMRNASeq[seqPos] == 'U') {
            pTotalScore += score;
        }

        pMaxScore += score;
    }

}

//
// TS5FeatThreePrimePair methods
//
void TS5FeatThreePrimePair::clear_features() {
    clear(mThreePrimePair);
    mAlign.clear_alignments();
}

int TS5FeatThreePrimePair::add_features(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        mikan::TSitePosSet const &pMRNAPos,
        mikan::TSitePosSet const &pSitePos,
        TS5FeatSeedType &pSeedTypes) {
    resize(mThreePrimePair, length(pMRNAPos));
    mAlign.resize_alignments((unsigned) length(pMRNAPos));
    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            mThreePrimePair[i] = -1.0;
            continue;
        }

        // Alias
        CharString const &seedType = pSeedTypes.get_seed_type(i);

        mThreePrimePair[i] = findBestMatch(i, pMRNAPos, pSitePos, seedType, pMRNASeqs[pMRNAPos[i]], pMiRNASeq);
        mAlign.align_seed(i, seedType, pMiRNASeq, pMRNASeqs[pMRNAPos[i]], (unsigned) pSitePos[i]);
    }

    return 0;
}

void TS5FeatThreePrimePair::getMRNASeq(
        CharString const &pSeedType,
        mikan::TRNAStr const &pMRNASeq,
        int pSitePos,
        mikan::TRNAStr &pMRNAThreePrime) {
    int startUTR = 0;
    int seqLen = 15;

    clear(pMRNAThreePrime);

    if (pSeedType == "8mer" || pSeedType == "7mer-m8") {
        startUTR = pSitePos - 16;
    } else if (pSeedType == "7mer-A1") {
        startUTR = pSitePos - 15;
    }

    if (startUTR < 0) {
        seqLen = 15 + startUTR;
        startUTR = 0;
    }
    resize(pMRNAThreePrime, seqLen);

    for (int i = 0; i < seqLen; ++i) {
        pMRNAThreePrime[i] = pMRNASeq[startUTR + i];
    }

    reverse(pMRNAThreePrime);

}

void TS5FeatThreePrimePair::getMiRNASeq(
        const CharString &pSeedType,
        const mikan::TRNAStr &pMiRNASeq,
        mikan::TRNAStr &pMiRNAThreePrime) {
    int startMiRNA = 0;
    int lenMiRNA;

    clear(pMiRNAThreePrime);

    if (pSeedType == "8mer" || pSeedType == "7mer-m8") {
        startMiRNA = 8;
    } else if (pSeedType == "7mer-A1") {
        startMiRNA = 7;
    }

    lenMiRNA = length(pMiRNASeq) - startMiRNA;
    resize(pMiRNAThreePrime, lenMiRNA);

    for (int i = 0; i < lenMiRNA; ++i) {
        pMiRNAThreePrime[i] = pMiRNASeq[startMiRNA + i];
    }

    complement(pMiRNAThreePrime);

}

float TS5FeatThreePrimePair::findBestMatch(
        unsigned pPosIdx,
        mikan::TSitePosSet const &,
        mikan::TSitePosSet const &pSitePos,
        const CharString &pSeedType,
        const mikan::TRNAStr &pMRNASeq,
        const mikan::TRNAStr &pMiRNASeq) {
    typedef Index<mikan::TRNAStr, IndexQGram<UngappedShape<2> > > TIndexQGram;
    typedef Finder<TIndexQGram> TFinder;
    mikan::TRNAStr mRNAThreePrime;
    mikan::TRNAStr miRNAThreePrime;
    // Get mRNA and miRNA seqs
    getMRNASeq(pSeedType, pMRNASeq, pSitePos[pPosIdx], mRNAThreePrime);
    getMiRNASeq(pSeedType, pMiRNASeq, miRNAThreePrime);

    float score = 0.0f;
    TIndexQGram RNAIdx(mRNAThreePrime);
    TFinder finder(RNAIdx);
    mikan::TRNAStr twoRNAs;
    String<int> matchLen;
    String<int> miRNAPos;
    String<int> mRNAPos;

    resize(twoRNAs, 2);

    for (unsigned i = 0; i < length(miRNAThreePrime) - 1; ++i) {
        twoRNAs[0] = miRNAThreePrime[i];
        twoRNAs[1] = miRNAThreePrime[i + 1];

        while (find(finder, twoRNAs)) {
            appendValue(matchLen, 2);
            appendValue(miRNAPos, i);
            appendValue(mRNAPos, position(finder));
        }

        goBegin(finder);
        clear(finder);
    }

    if (length(matchLen) > 1) {
        connectMatchedSeq(matchLen, miRNAPos, mRNAPos);
    }

    mIdxBestScore = 0;
    score = calcScore(pSeedType, matchLen, miRNAPos, mRNAPos);
    mAlign.align_3p_part(pPosIdx, pSeedType, miRNAThreePrime, mRNAThreePrime, matchLen, miRNAPos, mRNAPos,
                         score, mIdxBestScore);

    return score;
}

void TS5FeatThreePrimePair::connectMatchedSeq(
        String<int> &pMatchLen,
        String<int> &pMiRNAPos,
        String<int> &pMRNAPos) {
    int k = 2;
    bool allConnected = false;

    while (!allConnected) {
        allConnected = true;

        for (unsigned i = 0; i < length(pMatchLen) - 1; ++i) {
            if (pMatchLen[i] != k) {
                continue;
            }
            for (unsigned j = i + 1; j < length(pMatchLen); ++j) {
                if ((pMatchLen[j] != k) || (pMiRNAPos[j] == pMiRNAPos[i])) {
                    continue;
                }
                if ((pMiRNAPos[j] == pMiRNAPos[i] + 1) && (pMRNAPos[j] > pMRNAPos[i] + 1)) {
                    break;
                }
                if ((pMiRNAPos[j] == pMiRNAPos[i] + 1) && (pMRNAPos[j] == pMRNAPos[i] + 1)) {
                    allConnected = false;
                    pMatchLen[i] = pMatchLen[i] + 1;
                    break;
                }
                if (pMiRNAPos[j] > pMiRNAPos[i] + 1) {
                    break;
                }
            }
        }
        ++k;
    }

}

float TS5FeatThreePrimePair::calcScore(
        const CharString &pSeedType,
        String<int> &pMatchLen,
        String<int> &pMiRNAPos,
        String<int> &pMRNAPos) {
    float score, bestScore;
    int offset;
    int miRNAPos;
    int minOffSet;

    bestScore = 0.0;
    minOffSet = (int) length(pMatchLen);
    for (unsigned i = 0; i < length(pMatchLen); ++i) {
        score = 0.0;
        for (int j = 0; j < pMatchLen[i]; ++j) {
            miRNAPos = pMiRNAPos[i] + j;

            if (pSeedType == "8mer" || pSeedType == "7mer-m8") {
                if (4 <= miRNAPos && miRNAPos <= 7) {
                    score += 1.0;
                } else {
                    score += 0.5;
                }
            } else if (pSeedType == "7mer-A1") {
                if (5 <= miRNAPos && miRNAPos <= 8) {
                    score += 1.0;
                } else {
                    score += 0.5;
                }
            }

        }

        offset = std::abs(pMiRNAPos[i] - pMRNAPos[i]);
        score = score - std::max(0.0f, ((offset - 2.0f) / 2.0f));
        if (bestScore < score) {
            bestScore = score;
            mIdxBestScore = i;
            minOffSet = offset;
        } else if (score != 0.0 && bestScore == score) {
            if (offset < minOffSet) {
                mIdxBestScore = i;
                minOffSet = offset;
            }
        }
    }

    return bestScore;
}

} // namespace ts5cs
