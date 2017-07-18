#include <seqan/seq_io.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "rh2_seed_site.hpp"     // RH2SeedSeqs, RH2SeedSites

using namespace seqan;

namespace rh2mfe {

//
// RH2SeedSeqs methods
//
void RH2SeedSeqs::set_flags(mikan::TCharSet &pSeedTypeDef) {
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

    if (pSeedTypeDef[0][0] == '6') {
        mNMerLab = "6mer";
    } else if (pSeedTypeDef[0][0] == '7') {
        mNMerLab = "7mer";
    }

    if (pSeedTypeDef[0][2] == 'G' && pSeedTypeDef[0][3] == 'U') {
        CharString GUT;
        CharString GUM;

        if (pSeedTypeDef[0][0] == '7') {
            mGUTLab = "7mer_GUT";
            mGUMLab = "7mer_GUM";
        } else if (pSeedTypeDef[0][0] == '6') {
            mGUTLab = "6mer_GUT";
            mGUMLab = "6mer_GUM";
        } else {
            mGUTLab = "GUT";
            mGUMLab = "GUM";
        }

        if (pSeedTypeDef[0] == "6mGU+" || pSeedTypeDef[0] == "7mGU+") {
            mSingleGU = true;
            mMultiGU = true;
            mMultiGUTLab = mGUTLab;
            mMultiGUMLab = mGUMLab;
        } else if (pSeedTypeDef[0] == "6mGU1" || pSeedTypeDef[0] == "7mGU1") {
            mSingleGU = true;
        }
    }
}

//
// RH2SeedSites methods
//
bool RH2SeedSites::set_new_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr &pMiRNASeq,
        mikan::TCharSet &pSeedTypeDef,
        seqan::CharString &pSeedType,
        int,
        bool pEffectiveSite) {
    
    bool matchM8, matchM1, gutM8, gutM1, gumM8, gumM1, isA8, isA1, noM8, noM1;
    CharString newSeedType = "";
    
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 8, noM8, matchM8, gutM8, gumM8, isA8);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 1, noM1, matchM1, gutM1, gumM1, isA1);

    if (pSeedTypeDef[0][0] == '6') {
        newSeedType = pSeedType;
    } else if (pSeedTypeDef[0][0] == '7') {
        if (matchM8) {
            newSeedType = pSeedType;
        } else if (pSeedType == "7mer" || pSeedTypeDef[0] == "7mGU+") {
            if (gutM8) {
                newSeedType = "7mer_GUT";
            } else if (gumM8) {
                newSeedType = "7mer_GUM";
            } else {
                newSeedType = pSeedType;
                pEffectiveSite = false;
            }
        } else {
            newSeedType = pSeedType;
            pEffectiveSite = false;
        }
    } else if (pSeedTypeDef[0] == "targetscan") {
        if (!noM1 && isA1 && matchM8) {
            newSeedType = "8mer";
        } else if (!noM1 && isA1) {
            newSeedType = "7mer-A1";
        } else if (matchM8) {
            newSeedType = "7mer-m8";
        } else {
            newSeedType = "6mer";
            pEffectiveSite = false;
        }
    }

    if (pEffectiveSite) {
        appendValue(mSeedTypes, newSeedType);
        appendValue(mMisMatchPos, 0);
        appendValue(mEffectiveSites, pEffectiveSite);
    }

    return pEffectiveSite;

}

} // namespace rh2mfe
