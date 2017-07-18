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

    bool IsA1, MatchM8;
    mikan::TRNAStr revMiRNASeq, miRNAM8, mRNAM8, miRNAM8c, mRNAA1;
    CharString newSeedType = "";
    bool noA1;

    revMiRNASeq = pMiRNASeq;
    miRNAM8 = revMiRNASeq[7];
    complement(revMiRNASeq);
    miRNAM8c = revMiRNASeq[7];

    mRNAM8 = mMRNASeqs[pMRNAPos][pSitePos - 1];

    if (pSitePos + mikan::SEEDLEN < length(mMRNASeqs[pMRNAPos])) {
        mRNAA1 = mMRNASeqs[pMRNAPos][pSitePos + mikan::SEEDLEN];
        noA1 = false;
    } else {
        mRNAA1 = 'A';
        noA1 = true;
    }

    if (pSeedTypeDef[0][0] == '6') {
        newSeedType = pSeedType;
    } else if (pSeedTypeDef[0][0] == '7') {
        if (miRNAM8c == mRNAM8) {
            newSeedType = pSeedType;
        } else if (pSeedType == "7mer" || pSeedTypeDef[0] == "7mGU+") {
            if ((miRNAM8 == 'G') && (mRNAM8 == 'U')) {
                newSeedType = "7mer_GUT";
            } else if ((miRNAM8 == 'U') && (mRNAM8 == 'G')) {
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
        if (miRNAM8 == mRNAM8) {
            MatchM8 = true;
        } else {
            MatchM8 = false;
        }

        if (!noA1 && mRNAA1 == 'A') {
            IsA1 = true;
        } else {
            IsA1 = false;
        }

        if (!noA1 && IsA1 && MatchM8) {
            newSeedType = "8mer";
        } else if (!noA1 && IsA1) {
            newSeedType = "7mer-A1";
        } else if (MatchM8) {
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
