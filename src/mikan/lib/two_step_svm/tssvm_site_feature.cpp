#include <math.h>
#include <tssvm_align.hpp>          // TSAlign
#include <tssvm_inst_template.hpp>  // TRNATYPE
#include <tssvm_site_feature.hpp>   // TSSVMRawFeatures, TSSVMFeatSeedType, TSSVMFeatSimilarity,

using namespace seqan;

namespace tssvm {

//
// TSSVMRawFeatures methods
//
template<class TRNAString>
void TSSVMRawFeatures<TRNAString>::clear_features() {
    clear(mEffectiveSites);
    mSeedTypes.clear_features();
    mSimilarities.clear_features();
    mAURichUp.clear_features();
    mAURichDown.clear_features();
    mSitePos.clear_features();
    mSeqMatch.clear_features();
    mA1Match.clear_features();
}

template<class TRNAString>
int TSSVMRawFeatures<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        TSAlign<TRNAString> const &pAlignSeqs,
        TRNAString const &pMiRNASeq,
        StringSet<TRNAString> const &pMRNASeqs) {
    resize(mEffectiveSites, length(pSeedSites.mEffectiveSites));
    for (unsigned i = 0; i < length(mEffectiveSites); ++i) {
        mEffectiveSites[i] = pSeedSites.mEffectiveSites[i];
    }

    mSeedTypes.add_features(pSeedSites, pAlignSeqs, pMiRNASeq, pMRNASeqs, mEffectiveSites);
    mSimilarities.add_features(pSeedSites, pAlignSeqs, pMiRNASeq, pMRNASeqs, mEffectiveSites);
    mAURichUp.add_features(pSeedSites, pAlignSeqs, pMiRNASeq, pMRNASeqs, mEffectiveSites);
    mAURichDown.add_features(pSeedSites, pAlignSeqs, pMiRNASeq, pMRNASeqs, mEffectiveSites);
    mSitePos.add_features(pSeedSites, pAlignSeqs, pMiRNASeq, pMRNASeqs, mEffectiveSites);
    mSeqMatch.add_features(pSeedSites, pAlignSeqs, pMiRNASeq, pMRNASeqs, mEffectiveSites);
    mA1Match.add_features(pSeedSites, pAlignSeqs, pMiRNASeq, pMRNASeqs, mEffectiveSites);

    return 0;
}

//
// TSSVMFeatSeedType methods
//
template<class TRNAString>
void TSSVMFeatSeedType<TRNAString>::clear_features() {
    clear(mSeedTypes);
}

template<class TRNAString>
int TSSVMFeatSeedType<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        TSAlign<TRNAString> const &,
        TRNAString const &,
        StringSet<TRNAString> const &,
        String<bool> &pEffectiveSites) {
    const StringSet<CharString> &seedTypes = pSeedSites.get_seed_types();

    resize(mSeedTypes, length(pEffectiveSites));
    for (unsigned i = 0; i < length(pEffectiveSites); ++i) {
        resize(mSeedTypes[i], 9, 0);

        if (!pEffectiveSites[i]) {
            continue;
        }

        if (seedTypes[i] == "8mer") {
            mSeedTypes[i][0] = 1;
        } else if (seedTypes[i] == "7mer-m8") {
            mSeedTypes[i][1] = 1;
        } else if (seedTypes[i] == "7mer-A1") {
            mSeedTypes[i][2] = 1;
        } else if (seedTypes[i] == "6mer") {
            mSeedTypes[i][3] = 1;
        } else if (seedTypes[i] == "GUM") {
            mSeedTypes[i][4] = 1;
        } else if (seedTypes[i] == "GUT") {
            mSeedTypes[i][5] = 1;
        } else if (seedTypes[i] == "LP") {
            mSeedTypes[i][6] = 1;
        } else if (seedTypes[i] == "BT") {
            mSeedTypes[i][7] = 1;
        } else if (seedTypes[i] == "BM") {
            mSeedTypes[i][8] = 1;
        }

    }

    return 0;
}

//
// TSSVMFeatSimilarity methods
//
template<class TRNAString>
void TSSVMFeatSimilarity<TRNAString>::clear_features() {
    clear(mSimilarities);
}

template<class TRNAString>
int TSSVMFeatSimilarity<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &,
        TSAlign<TRNAString> const &pAlignSeqs,
        TRNAString const &,
        StringSet<TRNAString> const &,
        String<bool> &pEffectiveSites) {
    StringSet<CharString> const &alignBars = pAlignSeqs.get_align_bars();
    float simAll, simSeed, sim3p, simAddtional, tmpScore;
    unsigned simAllCount;
    char barChar;

    resize(mSimilarities, length(pEffectiveSites));

    for (unsigned i = 0; i < length(pEffectiveSites); ++i) {
        resize(mSimilarities[i], 4, 0);

        if (!pEffectiveSites[i]) {
            continue;
        }

        simAllCount = length(alignBars[i]);
        simAll = 0;
        simSeed = 0;
        sim3p = 0;
        simAddtional = 0;
        for (unsigned j = 0; j < simAllCount; ++j) {
            barChar = alignBars[i][simAllCount - j - 1];

            if (barChar == '|') {
                tmpScore = 1.0f;
            } else if (barChar == ':') {
                tmpScore = 0.4f;
            } else {
                tmpScore = 0.0f;
            }

            simAll += tmpScore;
            if (j < 9) {
                simSeed += tmpScore;
            } else {
                sim3p += tmpScore;
            }
            if (j > 11 && j < 16) {
                simAddtional += tmpScore;
            }
        }

        mSimilarities[i][0] = roundf(simAll / simAllCount * 10000) / 10000;
        mSimilarities[i][1] = roundf(simSeed / 9 * 10000) / 10000;
        mSimilarities[i][2] = roundf(sim3p / (simAllCount - 9) * 10000) / 10000;
        mSimilarities[i][3] = simAddtional / 4;
    }

    return 0;
}

//
// TSSVMFeatAURichUp methods
//
template<class TRNAString>
void TSSVMFeatAURichUp<TRNAString>::clear_features() {
    clear(mAURichUp);
}

template<class TRNAString>
int TSSVMFeatAURichUp<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        TSAlign<TRNAString> const &,
        TRNAString const &,
        StringSet<TRNAString> const &pMRNASeqs,
        String<bool> &pEffectiveSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &siteS8Pos = pSeedSites.get_site_pos_s8();
    unsigned startPos, endPos;

    resize(mAURichUp, length(pEffectiveSites));
    for (unsigned i = 0; i < length(mAURichUp); ++i) {
        resize(mAURichUp[i], 30, 0);
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start and end positions for upstream and downstream
        getUpStreamPos((unsigned) siteS8Pos[i], startPos, endPos);

        int k = 0;
        for (unsigned j = startPos; j < endPos; ++j) {
            if (pMRNASeqs[mRNAPos[i]][j] == 'A' || pMRNASeqs[mRNAPos[i]][j] == 'U') {
                mAURichUp[i][30 - (endPos - startPos) + k] = 1;
            }
            ++k;
        }
    }

    return 0;
}

template<class TRNAString>
void TSSVMFeatAURichUp<TRNAString>::getUpStreamPos(
        unsigned pS8Pos,
        unsigned &pStartPos,
        unsigned &pEndPos) {
    int ePos;

    pEndPos = pS8Pos;

    ePos = pEndPos - 30;
    if (ePos < 0) {
        pStartPos = 0;
    } else {
        pStartPos = pEndPos - 30;
    }
}

//
// TSSVMFeatAURichDown methods
//
template<class TRNAString>
void TSSVMFeatAURichDown<TRNAString>::clear_features() {
    clear(mAURichDown);
}

template<class TRNAString>
int TSSVMFeatAURichDown<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        TSAlign<TRNAString> const &,
        TRNAString const &,
        StringSet<TRNAString> const &pMRNASeqs,
        String<bool> &pEffectiveSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &siteS1Pos = pSeedSites.get_site_pos_s1();
    unsigned startPos, endPos;

    resize(mAURichDown, length(pEffectiveSites));
    startPos = 0;
    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        resize(mAURichDown[i], 30, 0);
        if (!pEffectiveSites[i]) {
            continue;
        }

        // Get start and end positions for upstream and downstream
        getDownStreamPos(siteS1Pos[i], length(pMRNASeqs[mRNAPos[i]]), startPos, endPos);

        int k = 0;
        for (unsigned j = startPos; j < endPos; ++j) {
            if (pMRNASeqs[mRNAPos[i]][j] == 'A' || pMRNASeqs[mRNAPos[i]][j] == 'U') {
                mAURichDown[i][k] = 1;
            }
            ++k;
        }
    }

    return 0;
}

template<class TRNAString>
void TSSVMFeatAURichDown<TRNAString>::getDownStreamPos(
        unsigned pS1Pos,
        unsigned pSeqLen,
        unsigned &pStartPos,
        unsigned &pEndPos) {
    pStartPos = pS1Pos;
    if (pS1Pos < pSeqLen) {
        pStartPos = pS1Pos;
        if ((pStartPos + 30) < pSeqLen) {
            pEndPos = pStartPos + 30;
        } else {
            pEndPos = pSeqLen;
        }
    } else {
        pStartPos = pSeqLen;
        pEndPos = pSeqLen;
    }

}

//
// TSSVMFeatSitePos methods
//
template<class TRNAString>
void TSSVMFeatSitePos<TRNAString>::clear_features() {
    clear(mSitePos);
}

template<class TRNAString>
int TSSVMFeatSitePos<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        TSAlign<TRNAString> const &,
        TRNAString const &,
        StringSet<TRNAString> const &pMRNASeqs,
        String<bool> &pEffectiveSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();
    const String<unsigned> &siteS8Pos = pSeedSites.get_site_pos_s8();
    float posScore;


    resize(mSitePos, length(pEffectiveSites));
    for (unsigned i = 0; i < length(mSitePos); ++i) {
        resize(mSitePos[i], 1, 0);
        if (!pEffectiveSites[i]) {
            continue;
        }

        posScore = (float) (siteS8Pos[i] + 1);
        posScore = roundf(posScore / length(pMRNASeqs[mRNAPos[i]]) * 10000) / 10000;

        mSitePos[i][0] = posScore;
    }

    return 0;
}

//
// TSSVMFeatSeqMatch methods
//
template<class TRNAString>
void TSSVMFeatSeqMatch<TRNAString>::clear_features() {
    clear(mSeqMatch);
}

template<class TRNAString>
int TSSVMFeatSeqMatch<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &,
        TSAlign<TRNAString> const &pAlignSeqs,
        TRNAString const &,
        StringSet<TRNAString> const &,
        String<bool> &pEffectiveSites) {
    StringSet<CharString> const &alignBars = pAlignSeqs.get_align_bars();
    char barChar;
    unsigned maxSeqLen;

    resize(mSeqMatch, length(pEffectiveSites));

    for (unsigned i = 0; i < length(mSeqMatch); ++i) {
        resize(mSeqMatch[i], 20, 0);

        if (!pEffectiveSites[i]) {
            continue;
        }

        maxSeqLen = length(alignBars[i]);
        if (maxSeqLen > 20) {
            maxSeqLen = 20;
        }
        for (unsigned j = 0; j < maxSeqLen; ++j) {
            barChar = alignBars[i][length(alignBars[i]) - j - 1];

            if (barChar == '|') {
                mSeqMatch[i][j] = 1.0f;
            } else if (barChar == ':') {
                mSeqMatch[i][j] = 0.4f;
            } else {
                mSeqMatch[i][j] = 0.0f;
            }
        }
    }

    return 0;
}

//
// TSSVMFeatA1Match methods
//
template<class TRNAString>
void TSSVMFeatA1Match<TRNAString>::clear_features() {
    clear(mA1Match);
}

template<class TRNAString>
int TSSVMFeatA1Match<TRNAString>::add_features(
        TSSVMSeedSites<TRNAString> &pSeedSites,
        TSAlign<TRNAString> const &,
        TRNAString const &,
        StringSet<TRNAString> const &,
        String<bool> &pEffectiveSites) {
    const StringSet<CharString> &seedTypes = pSeedSites.get_seed_types();

    resize(mA1Match, length(pEffectiveSites));
    for (unsigned i = 0; i < length(pEffectiveSites); ++i) {
        resize(mA1Match[i], 1, 0);

        if (!pEffectiveSites[i]) {
            continue;
        }

        if (seedTypes[i] == "7mer-m8" || seedTypes[i] == "6mer") {
            mA1Match[i][0] = 0;
        } else {
            mA1Match[i][0] = 1;
        }

    }

    return 0;
}

// Explicit template instantiation
template
class TSSVMFeatSeedType<TRNATYPE>;

template
class TSSVMFeatSimilarity<TRNATYPE>;

template
class TSSVMFeatAURichUp<TRNATYPE>;

template
class TSSVMFeatAURichDown<TRNATYPE>;

template
class TSSVMFeatSitePos<TRNATYPE>;

template
class TSSVMFeatSeqMatch<TRNATYPE>;

template
class TSSVMFeatA1Match<TRNATYPE>;

template
class TSSVMRawFeatures<TRNATYPE>;

} // namespace tssvm
