#include <seqan/seq_io.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_seed_site.hpp"      // TSSVMSequences, TSSVMSeedSeqs, TSSVMSeedSites, TSSVMSiteFilter

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

//
// TSSVMSiteFilter methods
//
float TSSVMSiteFilter::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &) {

    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    seqan::CharString seedType = seedTypes[pSitePos];
    float preced;

    if (seedType == "8mer" || seedType == "7mer-A1" || seedType == "7mer-m8") {
        preced = 0;
    } else if (seedType == "6mer") {
        preced = 1;
    } else if (seedType == "GUM") {
        preced = 2;
    } else if (seedType == "GUT") {
        preced = 3;
    } else if (seedType == "BT") {
        preced = 4;
    } else if (seedType == "BM") {
        preced = 5;
    } else if (seedType == "LP") {
        preced = 6;
    } else {
        preced = 7;
    }

    return preced;
}

void TSSVMSiteFilter::set_intervals(
        mikan::MKSeedSites &pSeedSites,
        unsigned pSiteIdx,
        unsigned &pStartSearch,
        unsigned &pEndSearch,
        unsigned &pStartAdd,
        unsigned &pEndAdd,
        bool &pSearchOverlap) {

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();
    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();

    pStartSearch = 0;
    pEndSearch = 0;
    pStartAdd = sitePos[pSiteIdx] - 1;
    pEndAdd = sitePos[pSiteIdx] + 7;
    pSearchOverlap = true;

    if (seedTypes[pSiteIdx] == "8mer"
        || seedTypes[pSiteIdx] == "7mer-A1"
        || seedTypes[pSiteIdx] == "7mer-m8"
        || seedTypes[pSiteIdx] == "6mer") {
        pStartSearch = pStartAdd;
        pEndSearch = pEndAdd;
        pSearchOverlap = false;
    } else if (seedTypes[pSiteIdx] == "GUM"
               || seedTypes[pSiteIdx] == "GUT"
               || seedTypes[pSiteIdx] == "LP") {
        pStartSearch = pStartAdd - 1;
        pEndSearch = pEndAdd + 1;
    } else if (seedTypes[pSiteIdx] == "BT") {
        pStartSearch = pStartAdd - 1;
        pEndSearch = pEndAdd + 2;
    } else if (seedTypes[pSiteIdx] == "BM") {
        pStartAdd = sitePos[pSiteIdx];
        pEndAdd = sitePos[pSiteIdx] + 8;
        pStartSearch = pStartAdd - 1;
        pEndSearch = pEndAdd;
    }
}

} // namespace tssvm
