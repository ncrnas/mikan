#include <seqan/seq_io.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mr3_seed_site.hpp"      // MR3SeedSeqs, MR3SeedSites

using namespace seqan;

namespace mr3as {

//
// MR3SeedSeqs methods
//
void MR3SeedSeqs::set_flags(mikan::TCharSet &pSeedTypeDef) {
    mSingleGU = false;
    mMultiGU = false;
    mMisMatch = false;
    mGUMisMatch = false;
    mBT = false;
    mBTM8 = false;
    mBM = false;
    mLP = false;
    mOther = false;
    mAddInReverse = false;

    if (pSeedTypeDef[2] == 'Y' || pSeedTypeDef[1] == 'Y') {
        if (pSeedTypeDef[3] == '1' || pSeedTypeDef[3] == '+') {
            mSingleGU = true;
        }
        if (pSeedTypeDef[3] == '+') {
            mMultiGU = true;
        }
        if (pSeedTypeDef[4] != "0:0") {
            mMisMatch = true;
        }
        if (pSeedTypeDef[4] != "0:0" && (pSeedTypeDef[3] == '1' || pSeedTypeDef[3] == '+')) {
            mGUMisMatch = true;
        }

        if (pSeedTypeDef[5] == "1") {
            mBT = true;
        }

    }
}

//
// MR3SeedSites methods
//
bool MR3SeedSites::check_position(unsigned pMRNAPos, unsigned pSitePos, seqan::CharString &pSeedType) {
    bool effectiveSite;
    unsigned endPos;

    effectiveSite = true;
    endPos = pSitePos + INDEXED_SEQ_LEN;
    if ((endPos < MIN_DIST_TO_CDS) || (endPos + MIN_DIST_UTR_END > length(mMRNASeqs[pMRNAPos]))) {
        effectiveSite = false;
    }

    return effectiveSite;
}

bool MR3SeedSites::set_new_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr &pMiRNASeq,
        mikan::TCharSet &pSeedTypeDef,
        seqan::CharString &pSeedType,
        int pMisMatchPos,
        bool pEffectiveSite) {

    bool matchM8, matchM9, gutM8, gutM9, gumM8, gumM9;;
    CharString newSeedType = "";

    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 8, matchM8, gutM8, gumM8);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 9, matchM9, gutM9, gumM9);

    set_stringent_seed_type(pSeedType, pSeedTypeDef, matchM8, matchM9, pMisMatchPos, newSeedType);

    if (newSeedType == "" && (pSeedTypeDef[2] == 'Y' || pSeedTypeDef[1] == 'Y')) {
        if (pSeedTypeDef[3] == '1' || pSeedTypeDef[3] == '+') {
            set_single_gu_seed_type(pSeedType, pSeedTypeDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                                    pMisMatchPos, newSeedType);
        }

        if (newSeedType == "" && pSeedTypeDef[3] == '+') {
            set_multiple_gu_seed_type(pSeedType, pSeedTypeDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                                      pMisMatchPos, newSeedType);
        }

        if (newSeedType == "" && pSeedTypeDef[4] != "0:0") {
            set_mismatch_seed_type(pSeedType, pSeedTypeDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                                   pMisMatchPos, newSeedType);
        }

        if (newSeedType == "" && pSeedTypeDef[4] != "0:0" && (pSeedTypeDef[3] == '1' || pSeedTypeDef[3] == '+')) {
            set_gu_mismatch_seed_type(pSeedType, pSeedTypeDef, -1, -2, matchM8, matchM9, gutM8, gutM9, gumM8, gumM9,
                                      pMisMatchPos, newSeedType);
        }

        if (pSeedType == "BT" && newSeedType == "" && pSeedTypeDef[5] == "1") {
            set_bt_seed_type(pMRNAPos, pSitePos, pMiRNASeq, pMisMatchPos, newSeedType);
        }

    }

    set_6mer_seed_type(pSeedType, pSeedTypeDef, matchM8, matchM9, pMisMatchPos, newSeedType);

    if (newSeedType != "") {
        appendValue(mEffectiveSites, true);
        pEffectiveSite = true;
    } else {
        pEffectiveSite = false;
    }

    return pEffectiveSite;

}

void MR3SeedSites::set_mx_matches(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr const &pMiRNA,
        int pMx,
        bool &pMatchMx,
        bool &pGutMx,
        bool &pGumMx) {
    mikan::TRNAStr cMiRNASeq, miRNAMx, mRNAMx, miRNAMxC;

    miRNAMx = pMiRNA[pMx - 1];
    cMiRNASeq = pMiRNA;
    complement(cMiRNASeq);
    miRNAMxC = cMiRNASeq[pMx - 1];

    pMatchMx = false;
    pGutMx = false;
    pGumMx = false;

    if ((int) pSitePos - pMx + 1 + (int) INDEXED_SEQ_LEN < 0) {
        return;
    }
    mRNAMx = mMRNASeqs[pMRNAPos][pSitePos - (pMx - 1 - INDEXED_SEQ_LEN)];

    pMatchMx = false;
    if (miRNAMxC == mRNAMx) {
        pMatchMx = true;
    }

    pGutMx = false;
    pGumMx = false;
    if (miRNAMx == 'G' && mRNAMx == 'U') {
        pGutMx = true;
    } else if (miRNAMx == 'U' && mRNAMx == 'G') {
        pGumMx = true;
    }

}

void MR3SeedSites::set_stringent_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        bool pMatchMx8,
        bool pMatchMx9,
        unsigned,
        CharString &pNewSeedType) {
    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer") {
        if (pSeedDef[2] == 'Y' && pMatchMx8 && pMatchMx9) {
            pNewSeedType = "8mer";
        } else if (pSeedDef[1] == 'Y' && pMatchMx8) {
            pNewSeedType = "7mer";
        }
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

void MR3SeedSites::set_6mer_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        bool,
        bool,
        unsigned,
        CharString &pNewSeedType) {
    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer") {
        if (pSeedDef[0] == 'Y') {
            pNewSeedType = "6mer";
        }
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

void MR3SeedSites::set_single_gu_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx8,
        bool pMatchMx9,
        bool pGutMx8,
        bool pGutMx9,
        bool pGumMx8,
        bool pGumMx9,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer") {
        if (pSeedDef[2] == 'Y' && ((pMatchMx8 && pGutMx9) || (pMatchMx9 && pGutMx8))) {
            pNewSeedType = "8mer_GUT";
            if (pMatchMx8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pSeedDef[2] == 'Y' && ((pMatchMx8 && pGumMx9) || (pMatchMx9 && pGumMx8))) {
            pNewSeedType = "8mer_GUM";
            if (pMatchMx8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pSeedDef[1] == 'Y' && pGutMx8) {
            pNewSeedType = "7mer_GUT";
            mm = pM1;
        } else if (pSeedDef[1] == 'Y' && pGumMx8) {
            pNewSeedType = "7mer_GUM";
            mm = pM1;
        }
    } else if (pCurSeedType == "GUT") {
        if (pSeedDef[2] == 'Y' && pMatchMx8 && pMatchMx9) {
            pNewSeedType = "8mer_GUT";
            mm = pMisMatchPos;
        } else if (pSeedDef[1] == 'Y' && pMatchMx8) {
            pNewSeedType = "7mer_GUT";
            mm = pMisMatchPos;
        }
    } else if (pCurSeedType == "GUM") {
        if (pSeedDef[2] == 'Y' && pMatchMx8 && pMatchMx9) {
            pNewSeedType = "8mer_GUM";
            mm = pMisMatchPos;
        } else if (pSeedDef[1] == 'Y' && pMatchMx8) {
            pNewSeedType = "7mer_GUM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx8, pMatchMx9, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

void MR3SeedSites::set_multiple_gu_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int,
        int,
        bool pMatchMx8,
        bool pMatchMx9,
        bool pGutMx8,
        bool pGutMx9,
        bool pGumMx8,
        bool pGumMx9,
        unsigned,
        CharString &pNewSeedType) {
    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer" && pSeedDef[2] == 'Y' && (pGumMx8 || pGutMx8) && (pGumMx9 || pGutMx9)) {
        pNewSeedType = "8mer_GU+";
    } else if (pCurSeedType == "GUT" || pCurSeedType == "GUM") {
        if (pSeedDef[2] == 'Y' && (pGumMx8 || pGutMx8 || pMatchMx8) && (pGumMx9 || pGutMx9 || pMatchMx9)) {
            pNewSeedType = "8mer_GU+";
        } else if (pSeedDef[1] == 'Y' && (pGumMx8 || pGutMx8)) {
            pNewSeedType = "7mer_GU+";
        }
    } else if (pCurSeedType == "GU+") {
        if (pSeedDef[2] == 'Y' && (pGumMx8 || pGutMx8 || pMatchMx8) && (pGumMx9 || pGutMx9 || pMatchMx9)) {
            pNewSeedType = "8mer_GU+";
        } else if (pSeedDef[1] == 'Y' && (pGumMx8 || pGutMx8 || pMatchMx8)) {
            pNewSeedType = "7mer_GU+";
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx8, pMatchMx9, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

void MR3SeedSites::set_mismatch_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx8,
        bool pMatchMx9,
        bool pGutMx8,
        bool pGutMx9,
        bool pGumMx8,
        bool pGumMx9,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {

    int mm = 0;

    if (pNewSeedType != "") {
        return;
    }

    if (pSeedDef[2] == 'Y') {
        if (pCurSeedType == "6mer" && ((pMatchMx8 && !pMatchMx9 && !pGutMx9 && !pGumMx9)
                                       || (!pMatchMx8 && !pGutMx8 && !pGumMx8 && pMatchMx9))) {
            pNewSeedType = "8mer_MM";
            if (pMatchMx8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pCurSeedType == "MM" && pMatchMx8 && pMatchMx9) {
            pNewSeedType = "8mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1") {
        if (pCurSeedType == "6mer" && !pMatchMx8 && !pGutMx8 && !pGumMx8) {
            pNewSeedType = "7mer_MM";
            mm = pM1;
        } else if (pCurSeedType == "MM" && pMatchMx8) {
            pNewSeedType = "7mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx8, pMatchMx9, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

void MR3SeedSites::set_gu_mismatch_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx8,
        bool pMatchMx9,
        bool pGutMx8,
        bool pGutMx9,
        bool pGumMx8,
        bool pGumMx9,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pSeedDef[2] == 'Y') {
        if (pCurSeedType == "6mer" && ((!pMatchMx8 && !pMatchMx9 && (pGutMx8 || pGumMx8) && !(pGutMx9 || pGumMx9))
                                       || (!pMatchMx8 && !pMatchMx9 && !(pGutMx8 || pGumMx8) &&
                                           (pGutMx9 || pGumMx9)))) {
            pNewSeedType = "8mer_MMGU";
            if (pGutMx8 || pGumMx8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") &&
                   ((pMatchMx8 && !pMatchMx9 && !pGutMx9 && !pGumMx9)
                    || (!pMatchMx8 && pMatchMx9 && !pGutMx8 && !pGumMx8))) {
            pNewSeedType = "8mer_MMGU";
            if (pMatchMx8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pCurSeedType == "MM" && ((pMatchMx8 && (pGutMx9 || pGumMx9))
                                            || ((pGutMx8 || pGumMx8) && pMatchMx9))) {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        } else if (pCurSeedType == "MMGU" && pMatchMx8 && pMatchMx9) {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1") {
        if (pCurSeedType == "MM" && (pGutMx9 || pGumMx9)) {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        } else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") && (!pMatchMx8 && !pGutMx8 && !pGumMx8)) {
            pNewSeedType = "7mer_MMGU";
            mm = pM1;
        } else if (pCurSeedType == "MMGU" && pMatchMx8) {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx8, pMatchMx9, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

void MR3SeedSites::set_bt_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr const &pMiRNA,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    mikan::TRNAStr cMiRNASeq, miRNAM2C, miRNAM8C, miRNAM9C, mRNAM2, mRNAM7, mRNAM8;

    if (pNewSeedType != "") {
        return;
    }

    cMiRNASeq = pMiRNA;
    complement(cMiRNASeq);
    miRNAM2C = cMiRNASeq[1];
    miRNAM8C = cMiRNASeq[7];
    miRNAM9C = cMiRNASeq[8];

    if ((pSitePos + INDEXED_SEQ_LEN) >= length(mMRNASeqs[pMRNAPos])) {
        return;
    }
    mRNAM2 = mMRNASeqs[pMRNAPos][pSitePos + INDEXED_SEQ_LEN];

    if (((int) pSitePos - (7 - (int) INDEXED_SEQ_LEN)) < 0) {
        return;
    }
    mRNAM7 = mMRNASeqs[pMRNAPos][pSitePos - (7 - INDEXED_SEQ_LEN)];

    if (miRNAM2C == mRNAM2 && miRNAM8C == mRNAM7) {
        pNewSeedType = "7mer_BT";
    }

    if (pNewSeedType == "7mer_BT" && ((int) pSitePos - 8 + (int) INDEXED_SEQ_LEN >= 0)) {
        mRNAM8 = mMRNASeqs[pMRNAPos][pSitePos - (8 - INDEXED_SEQ_LEN)];
        if (miRNAM9C == mRNAM8) {
            pNewSeedType = "8mer_BT";
        }
    }

    //TODO Check BT sites - seed pos 3-7 instead of pos 2-6?
//    std::cout << std::endl;
//    std::cout << pMiRNA << std::endl;
//    std::cout << mRNAM2 << ","  << mRNAM7 << "," <<  mRNAM8 << std::endl;
//    std::cout << miRNAM2C << ","  << miRNAM8C << "," <<  miRNAM9C << std::endl;

    if (FORCE_LAST_MATCH && pNewSeedType != "8mer_BT") {
        pNewSeedType = "";
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, pMisMatchPos);
    }

    return;
}

void MR3SeedSites::check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType) {
    if (pNewSeedType == "6mer" || pNewSeedType == "7mer" || pNewSeedType == "8mer") {
        return;
    }

    if ((pNewSeedType[0] == '7' && !pMatchM8) || (pNewSeedType[0] == '8' && !pMatchM9)) {
        pNewSeedType = "";
    }
}

} // namespace mr3as
