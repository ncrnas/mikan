#include <seqan/seq_io.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "pita_seed_site.hpp"     // PITASequences, PITASeedSeqs, PITASeedSites

using namespace seqan;

namespace ptddg {

//
// PITASeedSeqs methods
//
void PITASeedSeqs::set_flags(mikan::TCharSet &pSeedTypeDef) {
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
            mGUTLab = "GU+";
            mGUMLab = "GU+";
        }
        if (pSeedTypeDef[4] != "0:0") {
            mMisMatch = true;
        }
        if (pSeedTypeDef[4] != "0:0" && (pSeedTypeDef[3] == '1' || pSeedTypeDef[3] == '+')) {
            mGUMisMatch = true;
        }
    }
}

//
// PITASeedSites methods
//
bool PITASeedSites::check_position(unsigned pMRNAPos, unsigned pSitePos) {
    bool effectiveSite;
    unsigned endPos;

    effectiveSite = true;
    endPos = pSitePos + INDEXED_SEQ_LEN;
    if ((endPos < MIN_DIST_TO_CDS) || (endPos + MIN_DIST_UTR_END > length(mMRNASeqs[pMRNAPos]))) {
        effectiveSite = false;
    }

    return effectiveSite;
}

void PITASeedSites::set_new_seed_type(
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

    }

    set_6mer_seed_type(pSeedType, pSeedTypeDef, matchM8, matchM9, pMisMatchPos, newSeedType);

    if (newSeedType == "") {
        pEffectiveSite = false;
        appendValue(mSeedTypes, pSeedType);
        appendValue(mMisMatchPos, 0);
    }

    appendValue(mEffectiveSites, pEffectiveSite);

    return;

}

void PITASeedSites::set_mx_matches(
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

    mRNAMx = mMRNASeqs[pMRNAPos][pSitePos - (pMx - 1 - INDEXED_SEQ_LEN)];

    pMatchMx = false;
    if (miRNAMxC == mRNAMx) {
        pMatchMx = true;
    }

    pGutMx = false;
    pGumMx = false;
    if ((miRNAMx == 'G' && mRNAMx == 'U')) {
        pGutMx = true;
    } else if ((miRNAMx == 'U' && mRNAMx == 'G')) {
        pGumMx = true;
    }

}

void PITASeedSites::set_stringent_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        bool pMatchMx1,
        bool pMatchMx2,
        unsigned,
        CharString &pNewSeedType) {
    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer") {
        if (pSeedDef[2] == 'Y' && pMatchMx1 && pMatchMx2) {
            pNewSeedType = "8mer";
        } else if (pSeedDef[1] == 'Y' && pMatchMx1) {
            pNewSeedType = "7mer";
        }
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

void PITASeedSites::set_6mer_seed_type(
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

void PITASeedSites::set_single_gu_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer") {
        if (pSeedDef[2] == 'Y' && ((pMatchMx1 && pGutMx2) || (pMatchMx2 && pGutMx1))) {
            pNewSeedType = "8mer_GUT";
            if (pMatchMx1) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pSeedDef[2] == 'Y' && ((pMatchMx1 && pGumMx2) || (pMatchMx2 && pGumMx1))) {
            pNewSeedType = "8mer_GUM";
            if (pMatchMx1) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pSeedDef[1] == 'Y' && pGutMx1) {
            pNewSeedType = "7mer_GUT";
            mm = pM1;
        } else if (pSeedDef[1] == 'Y' && pGumMx1) {
            pNewSeedType = "7mer_GUM";
            mm = pM1;
        }
    } else if (pCurSeedType == "GUT") {
        if (pSeedDef[2] == 'Y' && pMatchMx1 && pMatchMx2) {
            pNewSeedType = "8mer_GUT";
            mm = pMisMatchPos;
        } else if (pSeedDef[1] == 'Y' && pMatchMx1) {
            pNewSeedType = "7mer_GUT";
            mm = pMisMatchPos;
        }
    } else if (pCurSeedType == "GUM") {
        if (pSeedDef[2] == 'Y' && pMatchMx1 && pMatchMx2) {
            pNewSeedType = "8mer_GUM";
            mm = pMisMatchPos;
        } else if (pSeedDef[1] == 'Y' && pMatchMx1) {
            pNewSeedType = "7mer_GUM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

void PITASeedSites::set_multiple_gu_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int,
        int,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned,
        CharString &pNewSeedType) {
    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer" && pSeedDef[2] == 'Y' && (pGumMx1 || pGutMx1) && (pGumMx2 || pGutMx2)) {
        pNewSeedType = "8mer_GU+";
    } else if (pCurSeedType == "GUT" || pCurSeedType == "GUM") {
        if (pSeedDef[2] == 'Y' && (pGumMx1 || pGutMx1 || pMatchMx1) && (pGumMx2 || pGutMx2 || pMatchMx2)) {
            pNewSeedType = "8mer_GU+";
        } else if (pSeedDef[1] == 'Y' && (pGumMx1 || pGutMx1)) //TODO: Check pGumMx1 || pGutMx1 || pMatchMx1
        {
            pNewSeedType = "7mer_GU+";
        }
    } else if (pCurSeedType == "GU+") {
        if (pSeedDef[2] == 'Y' && (pGumMx1 || pGutMx1 || pMatchMx1) && (pGumMx2 || pGutMx2 || pMatchMx2)) {
            pNewSeedType = "8mer_GU+";
        } else if (pSeedDef[1] == 'Y' && (pGumMx1 || pGutMx1 || pMatchMx1)) {
            pNewSeedType = "7mer_GU+";
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, 0);
    }

    return;
}

void PITASeedSites::set_mismatch_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pSeedDef[2] == 'Y') {
        if (pCurSeedType == "6mer" && ((pMatchMx1 && !pMatchMx2 && !pGutMx2 && !pGumMx2)
                                       || (!pMatchMx1 && !pGutMx1 && !pGumMx1 && pMatchMx2))) {
            pNewSeedType = "8mer_MM";
            if (pMatchMx1) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pCurSeedType == "MM" && pMatchMx1 && pMatchMx2) {
            pNewSeedType = "8mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1") {
        if (pCurSeedType == "6mer" && !pMatchMx1 && !pGutMx1 &&
            !pGumMx1) //TODO: Check !pMatchMx1 && !pGutMx1 && !pGumMx1
        {
            pNewSeedType = "7mer_MM";
            mm = pM1;
        } else if (pCurSeedType == "MM" && pMatchMx1) {
            pNewSeedType = "7mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

void PITASeedSites::set_gu_mismatch_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        int pM1,
        int pM2,
        bool pMatchMx1,
        bool pMatchMx2,
        bool pGutMx1,
        bool pGutMx2,
        bool pGumMx1,
        bool pGumMx2,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pSeedDef[2] == 'Y') {
        if (pCurSeedType == "6mer" && ((!pMatchMx1 && !pMatchMx2 && (pGutMx1 || pGumMx1) && !(pGutMx2 || pGumMx2))
                                       || (!pMatchMx1 && !pMatchMx2 && !(pGutMx1 || pGumMx1) &&
                                           (pGutMx2 || pGumMx2)))) {
            pNewSeedType = "8mer_MMGU";
            if (pGutMx1 || pGumMx1) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") &&
                   ((pMatchMx1 && !pMatchMx2 && !pGutMx2 && !pGumMx2)
                    || (!pMatchMx1 && pMatchMx2 && !pGutMx1 && !pGumMx1))) {
            pNewSeedType = "8mer_MMGU";
            if (pMatchMx1) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pCurSeedType == "MM" && ((pMatchMx1 && (pGutMx2 || pGumMx2))
                                            || ((pGutMx1 || pGumMx1) && pMatchMx2))) {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        } else if (pCurSeedType == "MMGU" && pMatchMx1 && pMatchMx2) {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1") {
        if (pCurSeedType == "MM" && (pGutMx2 || pGumMx2)) {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        } else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") && (!pMatchMx1 && !pGutMx1 && !pGumMx1)) {
            pNewSeedType = "7mer_MMGU";
            mm = pM1;
        } else if (pCurSeedType == "MMGU" && pMatchMx1) {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        }

    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchMx1, pMatchMx2, pNewSeedType);
    }

    if (pNewSeedType != "") {
        appendValue(mSeedTypes, pNewSeedType);
        appendValue(mMisMatchPos, mm);
    }

    return;
}

void PITASeedSites::check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType) {
    if (pNewSeedType == "6mer" || pNewSeedType == "7mer" || pNewSeedType == "8mer") {
        return;
    }

    if ((pNewSeedType[0] == '7' && !pMatchM8) || (pNewSeedType[0] == '8' && !pMatchM9)) {
        pNewSeedType = "";
    }
}

} // namespace ptddg
