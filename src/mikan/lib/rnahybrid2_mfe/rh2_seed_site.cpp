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
void RH2SeedSites::reset_finder() {
    goBegin(mFinder);
    clear(mFinder);
}

int RH2SeedSites::find_seed_sites(mikan::TRNAStr const &pMiRNA, mikan::TCharSet &pSeedDef) {
    RH2SeedSeqs seedSeqs;
    mikan::TRNAStr seedSeq;
    CharString seedType;
    int retVal;
    unsigned mRNAPos, sitePos;
    bool effectiveSite;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    seedSeqs.set_flags(pSeedDef);
    retVal = seedSeqs.create_seed_seqs();
    if (retVal != 0) {
        std::cerr << "ERROR: Could not get the seed sequence for " << pMiRNA;
        std::cerr << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(seedSeqs.mEffectiveSeeds); ++i) {
        if (!seedSeqs.mEffectiveSeeds[i]) {
            continue;
        }

        seedSeq = seedSeqs.get_seed_seq(i);
        seedType = seedSeqs.get_seed_type(i);

        while (find(mFinder, seedSeq)) {
            mRNAPos = position(mFinder).i1;
            sitePos = position(mFinder).i2;

            effectiveSite = true;
            if ((sitePos < MIN_DIST_TO_CDS) || (sitePos + MIN_DIST_UTR_END > length(mMRNASeqs[mRNAPos]))) {
                effectiveSite = false;
            }

            appendValue(mMRNAPos, mRNAPos);
            appendValue(mSitePos, sitePos);
            set_new_seed_type(seedType, pSeedDef[0], mRNAPos, sitePos, pMiRNA, effectiveSite);
            appendValue(mEffectiveSites, effectiveSite);
        }
        reset_finder();
    }

    return 0;
}

void RH2SeedSites::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mEffectiveSites);
}

void RH2SeedSites::set_new_seed_type(
        CharString &pCurSeedType,
        CharString &pSeedDef,
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr const &pMiRNA,
        bool &pEffectiveSite) {
    bool IsA1, MatchM8;
    mikan::TRNAStr revMiRNASeq, miRNAM8, mRNAM8, miRNAM8c, mRNAA1;
    CharString newSeedType = "";
    mikan::TRNAStr miRNASeq = pMiRNA;
    bool noA1;

    if (!pEffectiveSite) {
        appendValue(mSeedTypes, "");
        return;
    }
    revMiRNASeq = pMiRNA;
    miRNAM8 = revMiRNASeq[7];
    complement(revMiRNASeq);
    miRNAM8c = revMiRNASeq[7];

    mRNAM8 = mMRNASeqs[pMRNAPos][pSitePos - 1];

    if (pSitePos + SEED_LEN < length(mMRNASeqs[pMRNAPos])) {
        mRNAA1 = mMRNASeqs[pMRNAPos][pSitePos + SEED_LEN];
        noA1 = false;
    } else {
        mRNAA1 = 'A';
        noA1 = true;
    }

    if (pSeedDef[0] == '6') {
        newSeedType = pCurSeedType;
    } else if (pSeedDef[0] == '7') {
        if (miRNAM8c == mRNAM8) {
            newSeedType = pCurSeedType;
        } else if (pCurSeedType == "7mer" || pSeedDef == "7mGU+") {
            if ((miRNAM8 == 'G') && (mRNAM8 == 'U')) {
                newSeedType = "7mer_GUT";
            } else if ((miRNAM8 == 'U') && (mRNAM8 == 'G')) {
                newSeedType = "7mer_GUM";
            } else {
                newSeedType = pCurSeedType;
                pEffectiveSite = false;
            }
        } else {
            newSeedType = pCurSeedType;
            pEffectiveSite = false;
        }
    } else if (pSeedDef == "targetscan") {
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

    appendValue(mSeedTypes, newSeedType);

    return;

}

} // namespace rh2mfe
