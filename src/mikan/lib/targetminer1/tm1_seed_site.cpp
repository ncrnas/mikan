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
void TM1SeedSites::reset_finder() {
    goBegin(mFinder);
    clear(mFinder);
}

int TM1SeedSites::find_seed_sites(mikan::TRNAStr const &pMiRNA) {
    TM1SeedSeqs seedSeqs;
    mikan::TRNAStr seedSeq;
    CharString seedType;
    int retVal;
    unsigned mRNAPos, sitePos;
    bool effectiveSite;
    unsigned idx;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    mikan::TCharSet mNullSet;
    seedSeqs.set_flags(mNullSet);
    retVal = seedSeqs.create_seed_seqs();
    if (retVal != 0) {
        std::cerr << "ERROR: Could not get the seed sequence for " << pMiRNA;
        std::cerr << std::endl;
        return 1;
    }

    idx = 0;
    for (unsigned i = 0; i < length(seedSeqs.mEffectiveSeeds); ++i) {
        if (!seedSeqs.mEffectiveSeeds[i]) {
            continue;
        }

        seedSeq = seedSeqs.get_seed_seq(i);
        seedType = seedSeqs.get_seed_type(i);

        while (find(mFinder, seedSeq)) {
            mRNAPos = position(mFinder).i1;
            sitePos = position(mFinder).i2;

            appendValue(mMRNAPos, mRNAPos);
            appendValue(mSitePos, sitePos);
            effectiveSite = false;
            set_new_seed_type(seedType, mRNAPos, sitePos, pMiRNA, effectiveSite);
            appendValue(mEffectiveSites, effectiveSite);
            ++idx;
        }
        reset_finder();
    }

    return 0;
}

void TM1SeedSites::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mEffectiveSites);
    clear(mM8Match);
    clear(mM8GU);
    clear(mM1A);
    clear(mM1Match);
    clear(mM1GU);
    clear(mMRNASeqLen);
    clear(mM8Pos);
}

void TM1SeedSites::get_mx_match(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNAStr const &pMiRNACompSeq,
        mikan::TRNAStr const &pMRNASeq,
        unsigned pSitePos,
        int pMx,
        bool &pMatch,
        bool &pGU,
        bool &pNoMx,
        bool &pIsA) {
    int pos;
    mikan::TRNAStr miRNAMx, miRNAMxC, mRNAMx;

    pos = pSitePos - (pMx - 7);
    miRNAMx = pMiRNASeq[pMx - 1];
    miRNAMxC = pMiRNACompSeq[pMx - 1];

    pMatch = false;
    pGU = false;
    pIsA = false;

    if (pos >= 0 && pos < (int) length(pMRNASeq)) {
        mRNAMx = pMRNASeq[pos];
        if (mRNAMx == 'A') {
            pIsA = true;
        }
        if (miRNAMxC == mRNAMx) {
            pMatch = true;
        }
        if (((miRNAMx == 'G') && (mRNAMx == 'U')) || ((miRNAMx == 'U') && (mRNAMx == 'G'))) {
            pGU = true;
        }
        pNoMx = false;
    } else {
        pNoMx = true;
    }
}

void TM1SeedSites::get_match_count(
        mikan::TRNAStr const &pMiRNASeq,
        mikan::TRNAStr const &pMiRNACompSeq,
        mikan::TRNAStr const &pMRNASeq,
        unsigned pSitePos,
        int pMx1,
        int pMx2,
        int &pMatchCount,
        int &pGUCount) {
    bool isAx, matchMx, guMx, noMx;
    pMatchCount = 0;
    pGUCount = 0;

    for (int i = pMx1; i < pMx2 + 1; ++i) {
        get_mx_match(pMiRNASeq, pMiRNACompSeq, pMRNASeq, pSitePos, i, matchMx, guMx, noMx, isAx);
        if (!noMx) {
            if (matchMx) {
                ++pMatchCount;
            } else if (guMx) {
                ++pGUCount;
            }
        }
    }
}

void TM1SeedSites::set_new_seed_type(
        CharString &,
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr const &pMiRNA,
        bool &pEffectiveSite) {
    CharString newSeedType = "";
    bool isA1, isAx;
    bool matchM1, matchM2, matchM8, matchM9;
    bool guM1, guM2, guM8, guM9;
    bool noM1, noM2, noM8, noM9;
    int matchCount, guCount, matchCount2, guCount2;
    mikan::TRNAStr revMiRNASeq;

    revMiRNASeq = pMiRNA;
    complement(revMiRNASeq);

    // check m8, m2, m1 match
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 9, matchM9, guM9, noM9, isAx);
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 8, matchM8, guM8, noM8, isAx);
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 2, matchM2, guM2, noM2, isAx);
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 1, matchM1, guM1, noM1, isA1);

    get_match_count(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 2, 7, matchCount, guCount);
    get_match_count(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 2, 8, matchCount2, guCount2);

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
        pEffectiveSite = true;
    }

    appendValue(mSeedTypes, newSeedType);

    appendValue(mM8Match, matchM8);
    appendValue(mM8GU, guM8);
    appendValue(mM1A, isA1);
    appendValue(mM1Match, matchM1);
    appendValue(mM1GU, guM1);

    appendValue(mMRNASeqLen, length(mMRNASeqs[pMRNAPos]));
    appendValue(mM8Pos, pSitePos - 1);
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

int TM1SeedSites::get_seed_start_pos(int pIdx) {
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

    return std::max((int) mM8Pos[pIdx] + offset, 0);
}

int TM1SeedSites::get_seed_start_pos2(int pIdx) {
    int offset = 1;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match_gu(pIdx)) {
            offset = 0;
        } else {
            offset = 1;
        }
    }

    return std::max((int) mM8Pos[pIdx] + offset, 0);
}

int TM1SeedSites::get_length_to_cds(int pIdx) {
    int offset = 0;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 0;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match(pIdx)) {
            offset = 0;
        } else {
            offset = 1;
        }
    }

    return std::max((int) mM8Pos[pIdx] + offset, 0);
}

int TM1SeedSites::get_seed_end_pos(int pIdx) {
    int offset = 0;

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


    return std::min((int) mM8Pos[pIdx] + offset, (int) mMRNASeqLen[pIdx]);
}

int TM1SeedSites::get_seed_end_pos2(int pIdx) {
    int offset = 0;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 7;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match(pIdx)) {
            offset = 6;
        } else {
            offset = 7;
        }
    }


    return std::min((int) mM8Pos[pIdx] + offset, (int) mMRNASeqLen[pIdx]);
}

void TM1SeedSites::print_all() {
    for (unsigned i = 0 ; i < length(mEffectiveSites); i++) {
        std::cout << i << ", ";
        std::cout << mSeedTypes[i] << ", ";
        std::cout << mMRNAPos[i] << ", ";
        std::cout << mSitePos[i] << ", ";
        std::cout << mEffectiveSites[i] << std::endl;
    }

}

} // namespace tm1p
