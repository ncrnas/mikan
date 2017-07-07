#include <cmath>                 // min, max, round
#include "mk_typedef.hpp"        // TRNATYPE
#include "tm1_seed_site.hpp"     // TM1SeedSites
#include "tm1_site_feature.hpp"  // TM1RawFeatures, TM1FeatSeedType, TM1FeatSitePos, TM1FeatDistance, TM1FeatAURich,


using namespace seqan;
using namespace mikan;

namespace tm1p {

//
// TM1RawFeatures methods
//
template<class TRNAString>
void TM1RawFeatures<TRNAString>::clear_features() {
    mSeedTypes.clear_features();
    mSitePos.clear_features();
    mDistance.clear_features();
    mAURich.clear_features();
    mSingleFreqs.clear_features();
    mSingleFreqFlanks.clear_features();
    mDiFreqs.clear_features();
    mDiFreqFlanks.clear_features();
    mSingleMatches.clear_features();
    mTwoMatches.clear_features();

    mAlign.clear_alignments();
}

template<class TRNAString>
int TM1RawFeatures<TRNAString>::add_features(
        TRNAString const &pMiRNASeq,
        TRNASet const &pMRNASeqs,
        TM1SeedSites &pSeedSites,
        TM1SortedSitePos<TRNAString> &pSortedSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();

    mSeedTypes.add_features(pMiRNASeq, pMRNASeqs, pSeedSites.mEffectiveSites, pSeedSites);
    mSitePos.add_features(pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites, sitePos);
    mDistance.add_features(pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites, pSortedSites);
    mAURich.add_features(pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites, mDistance);
    mSingleFreqs.add_features(pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites);
    mSingleFreqFlanks.add_features(pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites, mDistance);
    mDiFreqs.add_features(pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites);
    mDiFreqFlanks.add_features(pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites, mDistance);
    mSingleMatches.add_features(pMiRNASeq, pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites, mSitePos);
    mTwoMatches.add_features(pMiRNASeq, pMRNASeqs, pSeedSites.mEffectiveSites, mRNAPos, pSeedSites, mSitePos);

    mAlign.align_all(pMiRNASeq, pMRNASeqs, pSeedSites, mSitePos.mM1Pos);

//    print_features(pSeedSites);

    return 0;
}

template<class TRNAString>
void TM1RawFeatures<TRNAString>::print_features(TM1SeedSites &pSeedSites) {
    for (unsigned i = 0; i < length(pSeedSites.mEffectiveSites); ++i) {
        if (!pSeedSites.mEffectiveSites[i]) {
            continue;
        }

        mSeedTypes.print_feature(i);
        mSitePos.print_feature(i);
        mDistance.print_feature(i);
        mAURich.print_feature(i);
        mSingleFreqs.print_feature(i);
        mSingleFreqFlanks.print_feature(i);
        mDiFreqs.print_feature(i);
        mDiFreqFlanks.print_feature(i);
        mSingleMatches.print_feature(i);
        mTwoMatches.print_feature(i);

        std::cout << std::endl;
    }

}

//
// TM1FeatSeedType methods
//
template<class TRNAString>
void TM1FeatSeedType<TRNAString>::clear_features() {
    clear(mSeedTypes);
}

template<class TRNAString>
void TM1FeatSeedType<TRNAString>::resize_features(unsigned pSize) {
    resize(mSeedTypes, pSize);
}

template<class TRNAString>
int TM1FeatSeedType<TRNAString>::add_features(
        TRNAString const &pMiRNASeq,
        TRNASet const &,
        String<bool> &pEffectiveSites,
        TM1SeedSites &pSeedSites) {
    TRNAString revMiRNASeq;
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();

    const StringSet<seqan::CharString> &seedTypes = pSeedSites.get_seed_types();

    resize_features((unsigned) length(mRNAPos));

    revMiRNASeq = pMiRNASeq;
    complement(revMiRNASeq);

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            mSeedTypes[i] = "";
            continue;
        }

        mSeedTypes[i] = seedTypes[i];
    }

    return 0;
}

template<class TRNAString>
void TM1FeatSeedType<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx << ", seed_type:" << mSeedTypes[pIdx] << std::endl;
}

//
// TM1FeatSitePos methods
//
template<class TRNAString>
void TM1FeatSitePos<TRNAString>::clear_features() {
    clear(mM1Pos);
    clear(mM8Pos);
    clear(mSeqLen);
}

template<class TRNAString>
void TM1FeatSitePos<TRNAString>::resize_features(unsigned pSize) {
    resize(mM1Pos, pSize);
    resize(mM8Pos, pSize);
    resize(mSeqLen, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        mM1Pos[i] = -1;
        mM8Pos[i] = -1;
        mSeqLen[i] = -1;
    }
}

template<class TRNAString>
int TM1FeatSitePos<TRNAString>::add_features(
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites,
        TSitePos const &pSitePos) {
    int seqLen;
    int relEndPosInt;
    float relEndPos;
    unsigned lenToCDS, endPos;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get m1 pos
        mM1Pos[i] = pSitePos[i] + SEARCH_SEQ_LEN;
        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        mSeqLen[i] = seqLen;

        // Check length to CDS
        mM8Pos[i] = pSitePos[i] - 1;
        lenToCDS = (unsigned) pSeedSites.get_seed_start_pos2(i);
        if (lenToCDS < MIN_DIST_TO_CDS || mM1Pos[i] > seqLen) {
            pEffectiveSites[i] = false;
            continue;
        }

        // Get relative pos
        endPos = (unsigned) pSeedSites.get_seed_end_pos2(i) - 1;
        relEndPos = (float) endPos / seqLen;
        relEndPos = roundf(relEndPos * 100.0f) / 100.0f;
        relEndPosInt = (int) (relEndPos * 1000);
        if (relEndPosInt >= 400 && relEndPosInt <= 600) {
            pEffectiveSites[i] = false;
            continue;
        }
    }

    return 0;
}

template<class TRNAString>
void TM1FeatSitePos<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx << ", site_pos:" << mM1Pos[pIdx] << ", seq_len:" << mSeqLen[pIdx] << std::endl;
}

//
// TM1FeatDistance methods
//
template<class TRNAString>
void TM1FeatDistance<TRNAString>::clear_features() {
    clear(mUpstream);
    clear(mDownstream);
}

template<class TRNAString>
void TM1FeatDistance<TRNAString>::resize_features(unsigned pSize) {
    resize(mUpstream, pSize);
    resize(mDownstream, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        mUpstream[i] = 0;
        mDownstream[i] = 0;
    }
}

template<class TRNAString>
int TM1FeatDistance<TRNAString>::add_features(
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites,
        TM1SortedSitePos<TRNAString> &pSortedSites) {
    const StringSet<String<unsigned> > &sortedSites = pSortedSites.get_sorted_mrna_pos();
    const String<unsigned> &mRNAIDs = pSortedSites.get_mrna_ids();
    int seedEnd, seedStart, seqLen, prevSeedEnd, prevSeedStart, idx;
    unsigned siteId;
    int startU, endU, lenU, startD, endD, lenD;
    String<int> nextSeedStarts;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(mRNAIDs); ++i) {
        seqLen = (int) length(pMRNASeqs[mRNAIDs[i]]);
        prevSeedStart = seqLen;
        resize(nextSeedStarts, length(sortedSites[i]));
        for (unsigned j = 0; j < length(sortedSites[i]); ++j) {
            idx = (int) length(sortedSites[i]) - j - 1;
            siteId = sortedSites[i][idx];
            if (!pEffectiveSites[siteId]) {
                nextSeedStarts[idx] = seqLen;
                continue;
            }
            nextSeedStarts[idx] = prevSeedStart;
            prevSeedStart = pSeedSites.get_seed_start_pos(siteId) + 1;
        }

        prevSeedEnd = 0;
        for (unsigned j = 0; j < length(sortedSites[i]); ++j) {
            siteId = sortedSites[i][j];
            if (!pEffectiveSites[siteId]) {
                continue;
            }

            seedStart = pSeedSites.get_seed_start_pos(siteId);
            seedEnd = pSeedSites.get_seed_end_pos(siteId);

            endU = pSeedSites.get_seed_start_pos(siteId);
            startU = std::max(0, endU - 30);
            startU = std::max(prevSeedEnd, startU);
            lenU = std::max(endU - startU, 0);
            mUpstream[siteId] = lenU;

            startD = pSeedSites.get_seed_end_pos(siteId);
            prevSeedEnd = startD - 1;
            endD = std::min(startD + 30, seqLen);
            endD = std::min(nextSeedStarts[j], endD);
            lenD = std::max(endD - startD, 0);
            mDownstream[siteId] = lenD;

        }
    }

    return 0;
}

template<class TRNAString>
void TM1FeatDistance<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx << ", upstream:" << mUpstream[pIdx];
    std::cout << ", downstream:" << mDownstream[pIdx] << std::endl;
}

//
// TM1FeatAURich methods
//
template<class TRNAString>
void TM1FeatAURich<TRNAString>::clear_features() {
    clear(mAURich);
}

template<class TRNAString>
void TM1FeatAURich<TRNAString>::resize_features(unsigned pSize) {
    resize(mAURich, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        mAURich[i] = -1.0;
    }
}

template<class TRNAString>
int TM1FeatAURich<TRNAString>::add_features(
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites,
        TM1FeatDistance<TRNAString> &pDistance) {
    (void) pDistance;
    int seqLen, startU, endU, startD, endD;
    float totalScore, upScore, upMaxScore, downScore, downMaxScore;
    //CharString chrUp = "up";
    //CharString chrDown = "down";

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start and end positions for upstream and downstream
        endU = pSeedSites.get_seed_start_pos(i);
        startU = std::max(0, endU - 30);

        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        startD = pSeedSites.get_seed_end_pos(i);
        endD = std::min(startD + 30, seqLen);

        calc_pos_scores(pMRNASeqs[pMRNAPos[i]], startU, endU, upScore, upMaxScore);
        calc_pos_scores(pMRNASeqs[pMRNAPos[i]], startD, endD, downScore, downMaxScore);

        totalScore = (upScore + downScore) / 60;

        mAURich[i] = totalScore;
    }

    return 0;
}

template<class TRNAString>
void TM1FeatAURich<TRNAString>::calc_pos_scores(
        const TRNAString &pMRNASeq,
        int pStart,
        int pEnd,
        float &pScore,
        float &pMaxScore) {
    int seqPos;

    pScore = 0;
    pMaxScore = 0;
    for (int i = 0; i < (pEnd - pStart); ++i) {
        seqPos = i + pStart;

        if (pMRNASeq[seqPos] == 'A' || pMRNASeq[seqPos] == 'U') {
            pScore += 1;
        }

        pMaxScore += 1;
    }
}

template<class TRNAString>
void TM1FeatAURich<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx << ", au_rich:" << mAURich[pIdx] << std::endl;
}

//
// TM1FeatSingleFreq methods
//
template<class TRNAString>
void TM1FeatSingleFreq<TRNAString>::clear_features() {
    clear(mSeedCountU);
    clear(mSeedCountC);
    clear(mSeedTotal);
}

template<class TRNAString>
void TM1FeatSingleFreq<TRNAString>::resize_features(unsigned pSize) {
    resize(mSeedCountU, pSize);
    resize(mSeedCountC, pSize);
    resize(mSeedTotal, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        mSeedCountU[i] = 0;
        mSeedCountC[i] = 0;
        mSeedTotal[i] = 0;
    }
}

template<class TRNAString>
int TM1FeatSingleFreq<TRNAString>::add_features(
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites) {
    int seedStart, seedEnd, curpos, seqLen;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start position for seed
        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        seedEnd = pSeedSites.get_seed_end_pos(i);
        seedStart = pSeedSites.get_seed_start_pos(i);

        for (int j = 0; j < (seedEnd - seedStart); ++j) {
            curpos = seedStart + j;

            if (pMRNASeqs[pMRNAPos[i]][curpos] == 'U') {
                ++mSeedCountU[i];
            } else if (pMRNASeqs[pMRNAPos[i]][curpos] == 'C') {
                ++mSeedCountC[i];
            }

        }

        mSeedTotal[i] = pSeedSites.get_seed_len(i);

    }


    return 0;
}

template<class TRNAString>
void TM1FeatSingleFreq<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx;
    std::cout << ", 14:" << mSeedCountU[pIdx];
    std::cout << ", 16:" << mSeedCountC[pIdx];
    std::cout << ", total_len:" << mSeedTotal[pIdx] << std::endl;
}

//
// TM1FeatSingleFreqFlank methods
//
template<class TRNAString>
void TM1FeatSingleFreqFlank<TRNAString>::clear_features() {
    clear(m3PCountA);
    clear(m3PCountU);
    clear(m3PTotal);
}

template<class TRNAString>
void TM1FeatSingleFreqFlank<TRNAString>::resize_features(unsigned pSize) {
    resize(m3PCountA, pSize);
    resize(m3PCountU, pSize);
    resize(m3PTotal, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        m3PCountA[i] = 0;
        m3PCountU[i] = 0;
        m3PTotal[i] = 0;
    }
}

template<class TRNAString>
int TM1FeatSingleFreqFlank<TRNAString>::add_features(
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites,
        TM1FeatDistance<TRNAString> &pDistance) {
    int seqLen, startU, endU, startD, endD;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {

        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start and end positions for upstream
        endU = pSeedSites.get_seed_start_pos(i);
        startU = std::max(0, endU - pDistance.get_upstream_len(i));
        calc_pos_scores(pMRNASeqs[pMRNAPos[i]], startU, endU, i);

        // Get start and end positions for downstream
        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        startD = pSeedSites.get_seed_end_pos(i);
        endD = std::min(startD + pDistance.get_downstream_len(i), seqLen);
        calc_pos_scores(pMRNASeqs[pMRNAPos[i]], startD, endD, i);

    }


    return 0;
}

template<class TRNAString>
void TM1FeatSingleFreqFlank<TRNAString>::calc_pos_scores(
        const TRNAString &pMRNASeq,
        int pSiteStart,
        int pSiteEnd,
        int pIdx) {
    int curpos;

    for (int i = 0; i < (pSiteEnd - pSiteStart); ++i) {
        curpos = pSiteStart + i;

        if (pMRNASeq[curpos] == 'A') {
            ++m3PCountA[pIdx];
        } else if (pMRNASeq[curpos] == 'U') {
            ++m3PCountU[pIdx];
        }
        ++m3PTotal[pIdx];
    }
}

template<class TRNAString>
void TM1FeatSingleFreqFlank<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx;
    std::cout << ", 17:" << m3PCountA[pIdx];
    std::cout << ", 18:" << m3PCountU[pIdx];
    std::cout << ", total_len:" << m3PTotal[pIdx] << std::endl;
}

//
// TM1FeatDiFreq methods
//
template<class TRNAString>
void TM1FeatDiFreq<TRNAString>::clear_features() {
    clear(mSeedCountUC);
    clear(mSeedCountCA);
    clear(mSeedCountCG);
    clear(mSeedTotal);
}

template<class TRNAString>
void TM1FeatDiFreq<TRNAString>::resize_features(unsigned pSize) {
    resize(mSeedCountUC, pSize);
    resize(mSeedCountCA, pSize);
    resize(mSeedCountCG, pSize);
    resize(mSeedTotal, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        mSeedCountUC[i] = 0;
        mSeedCountCA[i] = 0;
        mSeedCountCG[i] = 0;
        mSeedTotal[i] = 0;
    }
}

template<class TRNAString>
int TM1FeatDiFreq<TRNAString>::add_features(
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites) {
    int seedStart, seedEnd, seqLen, curpos1, curpos2;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start position for seed
        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        seedEnd = pSeedSites.get_seed_end_pos(i);
        seedStart = pSeedSites.get_seed_start_pos(i) - 1;

        for (int j = 0; j < (seedEnd - seedStart) - 1; ++j) {
            curpos1 = seedStart + j + 1;
            curpos2 = seedStart + j;

            if (pMRNASeqs[pMRNAPos[i]][curpos1] == 'U' && pMRNASeqs[pMRNAPos[i]][curpos2] == 'C') {
                ++mSeedCountUC[i];
            } else if (pMRNASeqs[pMRNAPos[i]][curpos1] == 'C' && pMRNASeqs[pMRNAPos[i]][curpos2] == 'A') {
                ++mSeedCountCA[i];
            } else if (pMRNASeqs[pMRNAPos[i]][curpos1] == 'C' && pMRNASeqs[pMRNAPos[i]][curpos2] == 'G') {
                ++mSeedCountCG[i];
            }
        }

        mSeedTotal[i] = pSeedSites.get_seed_len(i);

    }

    return 0;
}

template<class TRNAString>
void TM1FeatDiFreq<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx;
    std::cout << ", 28:" << mSeedCountUC[pIdx];
    std::cout << ", 33:" << mSeedCountCA[pIdx];
    std::cout << ", 35:" << mSeedCountCG[pIdx];
    std::cout << ", total_len:" << mSeedTotal[pIdx] << std::endl;
}

//
// TM1FeatDiFreqFlank methods
//
template<class TRNAString>
void TM1FeatDiFreqFlank<TRNAString>::clear_features() {
    clear(m3PCountAA);
    clear(m3PCountAU);
    clear(m3PCountAG);
    clear(m3PCountAC);
    clear(m3PCountUA);
    clear(m3PCountUU);
    clear(m3PCountUG);
    clear(m3PCountUC);
    clear(m3PCountGU);
    clear(m3PCountCA);
    clear(m3PCountCU);
    clear(m3PTotal);
}

template<class TRNAString>
void TM1FeatDiFreqFlank<TRNAString>::resize_features(unsigned pSize) {
    resize(m3PCountAA, pSize);
    resize(m3PCountAU, pSize);
    resize(m3PCountAG, pSize);
    resize(m3PCountAC, pSize);
    resize(m3PCountUA, pSize);
    resize(m3PCountUU, pSize);
    resize(m3PCountUG, pSize);
    resize(m3PCountUC, pSize);
    resize(m3PCountGU, pSize);
    resize(m3PCountCA, pSize);
    resize(m3PCountCU, pSize);
    resize(m3PTotal, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        m3PCountAA[i] = 0;
        m3PCountAU[i] = 0;
        m3PCountAG[i] = 0;
        m3PCountAC[i] = 0;
        m3PCountUA[i] = 0;
        m3PCountUU[i] = 0;
        m3PCountUG[i] = 0;
        m3PCountUC[i] = 0;
        m3PCountGU[i] = 0;
        m3PCountCA[i] = 0;
        m3PCountCU[i] = 0;
        m3PTotal[i] = 0;
    }
}

template<class TRNAString>
int TM1FeatDiFreqFlank<TRNAString>::add_features(
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites,
        TM1FeatDistance<TRNAString> &pDistance) {
    int seqLen, startU, endU, startD, endD, lenUp, lenDown;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start and end positions for upstream
        endU = pSeedSites.get_seed_start_pos(i);
        lenUp = pDistance.get_upstream_len(i);
        startU = std::max(0, endU - lenUp);
        calc_pos_scores(pMRNASeqs[pMRNAPos[i]], startU, endU, i);

        // Get start and end positions for downstream
        seqLen = (int) length(pMRNASeqs[pMRNAPos[i]]);
        startD = pSeedSites.get_seed_end_pos(i) - 1;
        lenDown = pDistance.get_downstream_len(i) + 1;
        endD = std::min(startD + lenDown, seqLen);
        calc_pos_scores(pMRNASeqs[pMRNAPos[i]], startD, endD, i);

        m3PTotal[i] = pDistance.get_upstream_len(i) + pDistance.get_downstream_len(i);
    }

    return 0;
}

template<class TRNAString>
void TM1FeatDiFreqFlank<TRNAString>::calc_pos_scores(
        const TRNAString &pMRNASeq,
        int pSiteStart,
        int pSiteEnd,
        int pIdx) {
    int curpos1, curpos2;

    for (int i = 0; i < (pSiteEnd - pSiteStart) - 1; ++i) {
        curpos1 = pSiteStart + i + 1;
        curpos2 = pSiteStart + i;

        if (pMRNASeq[curpos1] == 'A' && pMRNASeq[curpos2] == 'A') {
            ++m3PCountAA[pIdx];
        } else if (pMRNASeq[curpos1] == 'A' && pMRNASeq[curpos2] == 'U') {
            ++m3PCountAU[pIdx];
        } else if (pMRNASeq[curpos1] == 'A' && pMRNASeq[curpos2] == 'G') {
            ++m3PCountAG[pIdx];
        } else if (pMRNASeq[curpos1] == 'A' && pMRNASeq[curpos2] == 'C') {
            ++m3PCountAC[pIdx];
        } else if (pMRNASeq[curpos1] == 'U' && pMRNASeq[curpos2] == 'A') {
            ++m3PCountUA[pIdx];
        } else if (pMRNASeq[curpos1] == 'U' && pMRNASeq[curpos2] == 'U') {
            ++m3PCountUU[pIdx];
        } else if (pMRNASeq[curpos1] == 'U' && pMRNASeq[curpos2] == 'G') {
            ++m3PCountUG[pIdx];
        } else if (pMRNASeq[curpos1] == 'U' && pMRNASeq[curpos2] == 'C') {
            ++m3PCountUC[pIdx];
        } else if (pMRNASeq[curpos1] == 'G' && pMRNASeq[curpos2] == 'U') {
            ++m3PCountGU[pIdx];
        } else if (pMRNASeq[curpos1] == 'C' && pMRNASeq[curpos2] == 'A') {
            ++m3PCountCA[pIdx];
        } else if (pMRNASeq[curpos1] == 'C' && pMRNASeq[curpos2] == 'U') {
            ++m3PCountCU[pIdx];
        }

    }

}

template<class TRNAString>
void TM1FeatDiFreqFlank<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx;
    std::cout << ", 37:" << m3PCountAA[pIdx];
    std::cout << ", 38:" << m3PCountAU[pIdx];
    std::cout << ", 39:" << m3PCountAG[pIdx];
    std::cout << ", 40:" << m3PCountAC[pIdx];
    std::cout << ", 41:" << m3PCountUA[pIdx];
    std::cout << ", 42:" << m3PCountUU[pIdx];
    std::cout << ", 43:" << m3PCountUG[pIdx];
    std::cout << ", 44:" << m3PCountUC[pIdx];
    std::cout << ", 46:" << m3PCountGU[pIdx];
    std::cout << ", 49:" << m3PCountCA[pIdx];
    std::cout << ", 50:" << m3PCountCU[pIdx];
    std::cout << ", total_len:" << m3PTotal[pIdx] << std::endl;
}

//
// TM1FeatSingleMatch methods
//
template<class TRNAString>
void TM1FeatSingleMatch<TRNAString>::clear_features() {
    clear(mMatchUG);
    clear(mMatchGU);
    clear(mMatchCG);
    clear(mMatchTotal);
}

template<class TRNAString>
void TM1FeatSingleMatch<TRNAString>::resize_features(unsigned pSize) {
    resize(mMatchUG, pSize);
    resize(mMatchGU, pSize);
    resize(mMatchCG, pSize);
    resize(mMatchTotal, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        mMatchUG[i] = 0;
        mMatchGU[i] = 0;
        mMatchCG[i] = 0;
        mMatchTotal[i] = 0;
    }
}

template<class TRNAString>
int TM1FeatSingleMatch<TRNAString>::add_features(
        TRNAString const &pMiRNASeq,
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites,
        TM1FeatSitePos<TRNAString> &pSeedPos) {
    int seedStart, m1pos, curpos;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start position for seed
        m1pos = pSeedPos.get_m1_pos(i);
        seedStart = m1pos - 7;

        for (int j = 0; j < (m1pos - seedStart); ++j) {
            curpos = seedStart + j;

            if (pMRNASeqs[pMRNAPos[i]][curpos] == 'U' && pMiRNASeq[7 - j] == 'G') {
                ++mMatchUG[i];
                mMatchGU[i] = 0;
            } else if (pMRNASeqs[pMRNAPos[i]][curpos] == 'G' && pMiRNASeq[7 - j] == 'U') {
                mMatchUG[i] = 0;
                ++mMatchGU[i];
            } else if (pMRNASeqs[pMRNAPos[i]][curpos] == 'C' && pMiRNASeq[7 - j] == 'G') {
                ++mMatchCG[i];
            }

            if (mMatchUG[i] > 0) {
                mMatchUG[i] = 1;
            } else {
                mMatchUG[i] = 0;
            }

            if (mMatchGU[i] > 0) {
                mMatchGU[i] = 1;
            } else {
                mMatchGU[i] = 0;
            }

//            std::cout << i << "," << pMRNASeqs[pMRNAPos[i]][curpos] << "," << pMiRNASeq[7 - j] << std::endl;

        }

        mMatchTotal[i] = pSeedSites.get_seed_len(i);
    }


    return 0;
}

template<class TRNAString>
void TM1FeatSingleMatch<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx;
    std::cout << ", 55:" << mMatchUG[pIdx];
    std::cout << ", 57:" << mMatchGU[pIdx];
    std::cout << ", 58:" << mMatchCG[pIdx];
    std::cout << ", total_len:" << mMatchTotal[pIdx] << std::endl;
}

//
// TM1FeatTwoConsecMatch methods
//
template<class TRNAString>
void TM1FeatTwoConsecMatch<TRNAString>::clear_features() {
    clear(mMatchUACG);
    clear(mMatchUAUG);
    clear(mMatchCGGC);
    clear(mMatchUGGC);
    clear(mMatchTotal);
}

template<class TRNAString>
void TM1FeatTwoConsecMatch<TRNAString>::resize_features(unsigned pSize) {
    resize(mMatchUACG, pSize);
    resize(mMatchUAUG, pSize);
    resize(mMatchCGGC, pSize);
    resize(mMatchUGGC, pSize);
    resize(mMatchTotal, pSize);

    for (unsigned i = 0; i < pSize; ++i) {
        mMatchUACG[i] = 0;
        mMatchUAUG[i] = 0;
        mMatchCGGC[i] = 0;
        mMatchUGGC[i] = 0;
        mMatchTotal[i] = 1;
    }
}

template<class TRNAString>
int TM1FeatTwoConsecMatch<TRNAString>::add_features(
        TRNAString const &pMiRNASeq,
        TRNASet const &pMRNASeqs,
        String<bool> &pEffectiveSites,
        TSitePos const &pMRNAPos,
        TM1SeedSites &pSeedSites,
        TM1FeatSitePos<TRNAString> &pSeedPos) {
    int seedStart, m1pos, curpos1, curpos2;

    resize_features((unsigned) length(pMRNAPos));

    for (unsigned i = 0; i < length(pMRNAPos); ++i) {
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start position for seed
        m1pos = pSeedPos.get_m1_pos(i);
        seedStart = m1pos - 7;

        for (int j = 1; j < (m1pos - seedStart) - 1; ++j) {
            curpos1 = seedStart + j + 1;
            curpos2 = seedStart + j;

            if (pMRNASeqs[pMRNAPos[i]][curpos1] == 'U' && pMiRNASeq[6 - j] == 'A'
                && pMRNASeqs[pMRNAPos[i]][curpos2] == 'C' && pMiRNASeq[7 - j] == 'G') {
                ++mMatchUACG[i];
//                mMatchUAUG[i] = 0;
            } else if (pMRNASeqs[pMRNAPos[i]][curpos1] == 'U' && pMiRNASeq[6 - j] == 'A'
                       && pMRNASeqs[pMRNAPos[i]][curpos2] == 'U' && pMiRNASeq[7 - j] == 'G') {
//                mMatchUACG[i] = 0;
                ++mMatchUAUG[i];
            } else if (pMRNASeqs[pMRNAPos[i]][curpos1] == 'C' && pMiRNASeq[6 - j] == 'G'
                       && pMRNASeqs[pMRNAPos[i]][curpos2] == 'G' && pMiRNASeq[7 - j] == 'C') {
                ++mMatchCGGC[i];
            } else if (pMRNASeqs[pMRNAPos[i]][curpos1] == 'U' && pMiRNASeq[6 - j] == 'G'
                       && pMRNASeqs[pMRNAPos[i]][curpos2] == 'G' && pMiRNASeq[7 - j] == 'C') {
                ++mMatchUGGC[i];
            }

//            std::cout << i << "," << pMRNASeqs[pMRNAPos[i]][curpos1] <<  pMiRNASeq[6 - j];
//            std::cout << "-" << pMRNASeqs[pMRNAPos[i]][curpos2] << pMiRNASeq[7 - j] << std::endl;

            if (mMatchUAUG[i] > 0) {
                mMatchUAUG[i] = 1;
            } else {
                mMatchUAUG[i] = 0;
            }

            if (mMatchUAUG[i] > 0) {
                mMatchUAUG[i] = 1;
            } else {
                mMatchUAUG[i] = 0;
            }
        }

        mMatchTotal[i] = pSeedSites.get_seed_len(i);

    }


    return 0;
}

template<class TRNAString>
void TM1FeatTwoConsecMatch<TRNAString>::print_feature(unsigned pIdx) {
    std::cout << pIdx;
    std::cout << ", 68:" << mMatchUACG[pIdx];
    std::cout << ", 70:" << mMatchUAUG[pIdx];
    std::cout << ", 79:" << mMatchCGGC[pIdx];
    std::cout << ", 89:" << mMatchUGGC[pIdx];
    std::cout << ", total_len:" << mMatchTotal[pIdx] << std::endl;
}

// Explicit template instantiation
template
class TM1FeatSeedType<TRNATYPE>;

template
class TM1FeatSitePos<TRNATYPE>;

template
class TM1FeatDistance<TRNATYPE>;

template
class TM1FeatAURich<TRNATYPE>;

template
class TM1FeatSingleFreq<TRNATYPE>;

template
class TM1FeatSingleFreqFlank<TRNATYPE>;

template
class TM1FeatDiFreq<TRNATYPE>;

template
class TM1FeatDiFreqFlank<TRNATYPE>;

template
class TM1FeatSingleMatch<TRNATYPE>;

template
class TM1FeatTwoConsecMatch<TRNATYPE>;

template
class TM1RawFeatures<TRNATYPE>;

} // namespace tm1p
