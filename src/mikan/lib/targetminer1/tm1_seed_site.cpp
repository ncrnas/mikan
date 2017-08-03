#include <seqan/seq_io.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_seed_site.hpp"     // TM1SeedSeqs, TM1SeedSites

using namespace seqan;

namespace tm1p {

//
// TM1SeedSeqs methods
//
void TM1SeedSeqs::set_flags(mikan::TCharSet &) {
    mSingleGU = true;
    mMultiGU = false;
    mMisMatch = false;
    mGUMisMatch = false;
    mBT = false;
    mBTM8 = false;
    mBM = false;
    mLP = false;
    mOther = true;
    mAddInReverse = false;
}

int TM1SeedSeqs::create_other_seed_seqs(mikan::TRNAStr &pSeedSeq) {
    mikan::TRNAStr seedLPSeq;
    CharString seedType;
    unsigned m2pos = length(pSeedSeq) - 1;

    for (unsigned i = 0; i < length(mSeedSeqs); i++) {
        for (unsigned j = 0; j < length(mRNAChar); j++) {
            if (mSeedSeqs[i][m2pos] != mRNAChar[j]) {
                if (mSeedTypes[i] == "6mer") {
                    if ((mSeedSeqs[i][m2pos] == 'A' && mRNAChar[j] == 'G')
                        || (mSeedSeqs[i][m2pos] == 'C' && mRNAChar[j] == 'U')) {
                        continue;
                    }
                    seedType = "MM";
                } else if (mSeedTypes[i] == "GUM" || mSeedTypes[i] == "GUT") {
                    if ((mSeedSeqs[i][m2pos] == 'A' && mRNAChar[j] == 'G')
                        || (mSeedSeqs[i][m2pos] == 'C' && mRNAChar[j] == 'U')) {
                        seedType = "GUGU";
                    } else {
                        seedType = "GUMM";
                    }
                } else {
                    continue;
                }
                seedLPSeq = mSeedSeqs[i];
                seedLPSeq[m2pos] = mRNAChar[j];
                appendValue(mSeedSeqs, seedLPSeq);
                appendValue(mSeedTypes, seedType);
                appendValue(mMisMatchPos, m2pos);
            }
        }
    }

    return 0;
}

//
// TM1SeedSites methods
//
void TM1SeedSites::clear_pos() {
    mikan::MKSeedSites::clear_pos();

    clear(mM8Match);
    clear(mM8GU);
    clear(mM1A);
    clear(mM1Match);
    clear(mM1GU);
    clear(mMRNASeqLen);
}

void TM1SeedSites::get_match_count(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr const &pMiRNASeq,
        int pMx1,
        int pMx2,
        int &pMatchCount,
        int &pGUCount) {

    bool isAx, matchMx, gutMx, gumMx, noMx;
    pMatchCount = 0;
    pGUCount = 0;

    for (int i = pMx1; i < pMx2 + 1; ++i) {
        set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, i, noMx, matchMx, gutMx, gumMx, isAx);

        if (!noMx) {
            if (matchMx) {
                ++pMatchCount;
            } else if (gutMx || gumMx) {
                ++pGUCount;
            }
        }
    }
}

bool TM1SeedSites::check_position_2(unsigned pMRNAPos, unsigned pSitePos, seqan::CharString &pSeedType) {
    bool effectiveSite = true;
    unsigned curIdx = length(mEffectiveSites) - 1;
    unsigned m1Pos = pSitePos + INDEXED_SEQ_LEN;
    int m8Pos = pSitePos - 1;
    unsigned seqLen = length(mMRNASeqs[pMRNAPos]);
    unsigned offsetS = 1;
    unsigned offsetE = 7;
    unsigned lenToCDS, endPos;
    float relEndPos;
    int relEndPosInt;

    if (pSeedType == "7mer-m8") {
        offsetE = 6;
    } else if (pSeedType == "6mer") {
        if (is_m8_match_gu(curIdx)) {
            offsetS = 0;
        }
        if (is_m8_match(curIdx)) {
            offsetE = 5;
        } else {
            offsetE = 6;
        }
    }
    lenToCDS = (m8Pos + offsetS > 0) ? m8Pos + offsetS : 0;
    endPos = m8Pos + offsetE;


    if (lenToCDS < mMinToCDS || m1Pos > seqLen) {
        effectiveSite = false;
    }

    relEndPos = (float) endPos / seqLen;
    relEndPos = roundf(relEndPos * 100.0f) / 100.0f;
    relEndPosInt = (int) (relEndPos * 1000);
    if (relEndPosInt >= 400 && relEndPosInt <= 600) {
        effectiveSite = false;
    }

    return effectiveSite;
}

bool TM1SeedSites::set_new_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr &pMiRNASeq,
        mikan::TCharSet &,
        seqan::CharString &pSeedType,
        int,
        bool pEffectiveSite) {

    CharString newSeedType = "";
    bool isA1, isAx;
    bool matchM1, matchM2, matchM8, matchM9;
    bool gumM1, gutM2, gutM8, gutM9;
    bool gutM1, gumM2, gumM8, gumM9;
    bool noM1, noM2, noM8, noM9;
    int matchCount, guCount, matchCount2, guCount2;

    // check m8, m2, m1 match
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 9, noM9, matchM9, gutM9, gumM9, isAx);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 8, noM8, matchM8, gutM8, gumM8, isAx);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 2, noM2, matchM2, gutM2, gumM2, isAx);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 1, noM1, matchM1, gutM1, gumM1, isA1);

    get_match_count(pMRNAPos, pSitePos, pMiRNASeq, 2, 7, matchCount, guCount);
    get_match_count(pMRNAPos, pSitePos, pMiRNASeq, 2, 8, matchCount2, guCount2);

    if (isA1 && matchM8 && (matchCount + guCount == 6) && guCount < 2) {
        newSeedType = "8mer";
    } else if (!isA1 && matchM8 && (matchCount + guCount == 6) && guCount < 2) {
        newSeedType = "7mer-m8";
    } else if (isA1 && !matchM8 && (matchCount + guCount == 6) && guCount < 2) {
        newSeedType = "7mer-A1";
    } else if (((matchCount2 + guCount2 == 6) && guCount2 < 2) || ((matchCount2 + guCount2 > 6) && guCount2 < 3)) {
        newSeedType = "6mer";
    }

    if (newSeedType != "") {
        appendValue(mSeedTypes, newSeedType);
        pSeedType = newSeedType;

        appendValue(mM8Match, matchM8);
        appendValue(mM8GU, gutM8 || gumM8);
        appendValue(mM1A, isA1);
        appendValue(mM1Match, matchM1);
        appendValue(mM1GU, gutM1 || gumM1);

        appendValue(mMRNASeqLen, length(mMRNASeqs[pMRNAPos]));

        appendValue(mMisMatchPos, 0);

        appendValue(mEffectiveSites, true);

        pEffectiveSite = true;
    } else {
        pEffectiveSite = false;
    }

    return pEffectiveSite;

}

int TM1SeedSites::get_seed_len(int pIdx) {
    if (mSeedTypes[pIdx] == "8mer") {
        return 8;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        return 7;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        return 7;
    } else if (mSeedTypes[pIdx] == "6mer") {
        return 6;
    }

    return 6;
}

int TM1SeedSites::get_seed_start(int pIdx) {
    int offset = 0;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 0;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 0;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match_gu(pIdx)) {
            offset = 0;
        } else {
            offset = 1;
        }
    }

    return std::max((int) mS8Pos[pIdx] + offset, 0);

}

int TM1SeedSites::get_seed_end(int pIdx) {
    int offset = 0;
    int seqLen = get_mrna_seq_len(pIdx);

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 7;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match_gu(pIdx)) {
            offset = 6;
        } else {
            offset = 7;
        }

    }

    return std::min((int) mS8Pos[pIdx] + offset, seqLen);
}

} // namespace tm1p
