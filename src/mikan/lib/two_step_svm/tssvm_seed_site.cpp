#include <seqan/seq_io.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_seed_site.hpp"      // TSSVMSequences, TSSVMSeedSeqs, TSSVMSeedSites

using namespace seqan;

namespace tssvm {

//
// TSSVMSeedSeqs methods
//
void TSSVMSeedSeqs::set_flags(mikan::TCharSet &) {
    mSingleGU = true;
    mMultiGU = false;
    mMisMatch = false;
    mGUMisMatch = false;
    mBT = false;
    mBTM8 = true;
    mBM = true;
    mLP = true;
    mOther = true;
    mAddInReverse = true;
    mFilterRedundant = false;
    mTSSVMMismatch = true;

    mGUTLab = "GUM";
    mGUMLab = "GUT";
}

//
// TSSVMSeedSites methods
//
bool TSSVMSeedSites::set_new_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr &pMiRNASeq,
        mikan::TCharSet &,
        seqan::CharString &pSeedType,
        int pMisMatchPos,
        bool pEffectiveSite) {

    bool matchM8, matchM2, matchM1, gutMx, gumMx, isAx, isA1, noMx, noM1;
    CharString newSeedType = "";

    unsigned s8pos = pSitePos;
    unsigned s1pos = pSitePos;
    if (pSeedType == "BM") {
        s8pos = pSitePos + 1;
    } else if (pSeedType == "BT") {
        s1pos = pSitePos + 1;
    }
    set_mx_matches(pMRNAPos, s8pos, pMiRNASeq, 8, noMx, matchM8, gutMx, gumMx, isAx);
    set_mx_matches(pMRNAPos, s1pos, pMiRNASeq, 2, noMx, matchM2, gutMx, gumMx, isAx);
    set_mx_matches(pMRNAPos, s1pos, pMiRNASeq, 1, noM1, matchM1, gutMx, gumMx, isA1);

    if ((pSeedType != "6mer") && noM1) {
        return false;
    } else if (pSeedType == "6mer") {
        if ((isA1 || noM1) && matchM8) {
            newSeedType = "8mer";
        } else if (isA1 || noM1) {
            newSeedType = "7mer-A1";
        } else if (matchM8) {
            newSeedType = "7mer-m8";
        } else {
            newSeedType = "6mer";
        }
    } else if (isA1 && matchM8) {
        if (pSeedType == "BT") {
            if (matchM2) {
                newSeedType = pSeedType;
            }
        } else {
            newSeedType = pSeedType;
        }
    }

    if (newSeedType == "") {
        pEffectiveSite = false;
    } else {
        appendValue(mSeedTypes, newSeedType);
        appendValue(mMisMatchPos, pMisMatchPos);

        appendValue(mS1Pos, s1pos + 6);
        appendValue(mS8Pos, s8pos - 1);

        appendValue(mEffectiveSites, true);
        pEffectiveSite = true;
    }

    return pEffectiveSite;
}

} // namespace tssvm
