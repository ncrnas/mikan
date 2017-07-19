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
bool PITASeedSites::set_new_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr &pMiRNASeq,
        mikan::TCharSet &pSeedTypeDef,
        seqan::CharString &pSeedType,
        int pMisMatchPos,
        bool pEffectiveSite) {

    bool matchM8, matchM9, gutM8, gutM9, gumM8, gumM9, isAx, noMx;
    CharString newSeedType = "";

    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 8, noMx, matchM8, gutM8, gumM8, isAx);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 9, noMx, matchM9, gutM9, gumM9, isAx);

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

    if (newSeedType != "") {
        appendValue(mEffectiveSites, true);
        pEffectiveSite = true;
    } else {
        pEffectiveSite = false;
    }

    return pEffectiveSite;

}

void PITASeedSites::set_stringent_seed_type(
        CharString &pCurSeedType,
        StringSet<CharString> &pSeedDef,
        bool pMatchM8,
        bool pMatchM9,
        unsigned,
        CharString &pNewSeedType) {
    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer") {
        if (pSeedDef[2] == 'Y' && pMatchM8 && pMatchM9) {
            pNewSeedType = "8mer";
        } else if (pSeedDef[1] == 'Y' && pMatchM8) {
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
        bool pMatchM8,
        bool pMatchM9,
        bool pGutM8,
        bool pGutM9,
        bool pGumM8,
        bool pGumM9,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer") {
        if (pSeedDef[2] == 'Y' && ((pMatchM8 && pGutM9) || (pMatchM9 && pGutM8))) {
            pNewSeedType = "8mer_GUT";
            if (pMatchM8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pSeedDef[2] == 'Y' && ((pMatchM8 && pGumM9) || (pMatchM9 && pGumM8))) {
            pNewSeedType = "8mer_GUM";
            if (pMatchM8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pSeedDef[1] == 'Y' && pGutM8) {
            pNewSeedType = "7mer_GUT";
            mm = pM1;
        } else if (pSeedDef[1] == 'Y' && pGumM8) {
            pNewSeedType = "7mer_GUM";
            mm = pM1;
        }
    } else if (pCurSeedType == "GUT") {
        if (pSeedDef[2] == 'Y' && pMatchM8 && pMatchM9) {
            pNewSeedType = "8mer_GUT";
            mm = pMisMatchPos;
        } else if (pSeedDef[1] == 'Y' && pMatchM8) {
            pNewSeedType = "7mer_GUT";
            mm = pMisMatchPos;
        }
    } else if (pCurSeedType == "GUM") {
        if (pSeedDef[2] == 'Y' && pMatchM8 && pMatchM9) {
            pNewSeedType = "8mer_GUM";
            mm = pMisMatchPos;
        } else if (pSeedDef[1] == 'Y' && pMatchM8) {
            pNewSeedType = "7mer_GUM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchM8, pMatchM9, pNewSeedType);
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
        bool pMatchM8,
        bool pMatchM9,
        bool pGutM8,
        bool pGutM9,
        bool pGumM8,
        bool pGumM9,
        unsigned,
        CharString &pNewSeedType) {
    if (pNewSeedType != "") {
        return;
    }

    if (pCurSeedType == "6mer" && pSeedDef[2] == 'Y' && (pGumM8 || pGutM8) && (pGumM9 || pGutM9)) {
        pNewSeedType = "8mer_GU+";
    } else if (pCurSeedType == "GUT" || pCurSeedType == "GUM") {
        if (pSeedDef[2] == 'Y' && (pGumM8 || pGutM8 || pMatchM8) && (pGumM9 || pGutM9 || pMatchM9)) {
            pNewSeedType = "8mer_GU+";
        } else if (pSeedDef[1] == 'Y' && (pGumM8 || pGutM8)) //TODO: Check pGumM8 || pGutM8 || pMatchM8
        {
            pNewSeedType = "7mer_GU+";
        }
    } else if (pCurSeedType == "GU+") {
        if (pSeedDef[2] == 'Y' && (pGumM8 || pGutM8 || pMatchM8) && (pGumM9 || pGutM9 || pMatchM9)) {
            pNewSeedType = "8mer_GU+";
        } else if (pSeedDef[1] == 'Y' && (pGumM8 || pGutM8 || pMatchM8)) {
            pNewSeedType = "7mer_GU+";
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchM8, pMatchM9, pNewSeedType);
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
        bool pMatchM8,
        bool pMatchM9,
        bool pGutM8,
        bool pGutM9,
        bool pGumM8,
        bool pGumM9,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pSeedDef[2] == 'Y') {
        if (pCurSeedType == "6mer" && ((pMatchM8 && !pMatchM9 && !pGutM9 && !pGumM9)
                                       || (!pMatchM8 && !pGutM8 && !pGumM8 && pMatchM9))) {
            pNewSeedType = "8mer_MM";
            if (pMatchM8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pCurSeedType == "MM" && pMatchM8 && pMatchM9) {
            pNewSeedType = "8mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1") {
        if (pCurSeedType == "6mer" && !pMatchM8 && !pGutM8 &&
            !pGumM8) //TODO: Check !pMatchM8 && !pGutM8 && !pGumM8
        {
            pNewSeedType = "7mer_MM";
            mm = pM1;
        } else if (pCurSeedType == "MM" && pMatchM8) {
            pNewSeedType = "7mer_MM";
            mm = pMisMatchPos;
        }
    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchM8, pMatchM9, pNewSeedType);
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
        bool pMatchM8,
        bool pMatchM9,
        bool pGutM8,
        bool pGutM9,
        bool pGumM8,
        bool pGumM9,
        unsigned pMisMatchPos,
        CharString &pNewSeedType) {
    int mm;

    if (pNewSeedType != "") {
        return;
    }

    if (pSeedDef[2] == 'Y') {
        if (pCurSeedType == "6mer" && ((!pMatchM8 && !pMatchM9 && (pGutM8 || pGumM8) && !(pGutM9 || pGumM9))
                                       || (!pMatchM8 && !pMatchM9 && !(pGutM8 || pGumM8) &&
                                           (pGutM9 || pGumM9)))) {
            pNewSeedType = "8mer_MMGU";
            if (pGutM8 || pGumM8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") &&
                   ((pMatchM8 && !pMatchM9 && !pGutM9 && !pGumM9)
                    || (!pMatchM8 && pMatchM9 && !pGutM8 && !pGumM8))) {
            pNewSeedType = "8mer_MMGU";
            if (pMatchM8) {
                mm = pM2;
            } else {
                mm = pM1;
            }
        } else if (pCurSeedType == "MM" && ((pMatchM8 && (pGutM9 || pGumM9))
                                            || ((pGutM8 || pGumM8) && pMatchM9))) {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        } else if (pCurSeedType == "MMGU" && pMatchM8 && pMatchM9) {
            pNewSeedType = "8mer_MMGU";
            mm = pMisMatchPos;
        }
    }

    if (pNewSeedType == "" && pSeedDef[1] == 'Y' && pSeedDef[4] == "1:1") {
        if (pCurSeedType == "MM" && (pGutM9 || pGumM9)) {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        } else if ((pCurSeedType == "GUT" || pCurSeedType == "GUM") && (!pMatchM8 && !pGutM8 && !pGumM8)) {
            pNewSeedType = "7mer_MMGU";
            mm = pM1;
        } else if (pCurSeedType == "MMGU" && pMatchM8) {
            pNewSeedType = "7mer_MMGU";
            mm = pMisMatchPos;
        }

    }

    if (FORCE_LAST_MATCH && pNewSeedType != "") {
        check_last_match(pMatchM8, pMatchM9, pNewSeedType);
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
