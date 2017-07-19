#include <seqan/seq_io.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_seed_site.hpp"     // TS5SeedSeqs, TS5SeedSites

using namespace seqan;

namespace ts5cs {

//
// TS5SeedSeqs methods
//
void TS5SeedSeqs::set_flags(mikan::TCharSet &) {
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
}

//
// TS5SeedSites methods
//
bool TS5SeedSites::set_new_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr &pMiRNASeq,
        mikan::TCharSet &,
        seqan::CharString &,
        int pMisMatchPos,
        bool pEffectiveSite) {

    bool matchM8, matchM1, gutMx, gumMx, isAx, isA1, noMx, noM1;
    int lenToCds;
    CharString newSeedType = "";

    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 8, noMx, matchM8, gutMx, gumMx, isAx);
    set_mx_matches(pMRNAPos, pSitePos, pMiRNASeq, 1, noM1, matchM1, gutMx, gumMx, isA1);

    if (isA1 && matchM8) {
        newSeedType = "8mer";
    } else if (isA1) {
        newSeedType = "7mer-A1";
    } else if (matchM8) {
        newSeedType = "7mer-m8";
    }

    lenToCds = pSitePos;
    if (newSeedType == "8mer" || newSeedType == "7mer-m8") {
        lenToCds -= 1;
    }
    if (lenToCds < static_cast<int>(mMinToCDS)) {
        newSeedType = "";
    }

    if (newSeedType == "") {
        pEffectiveSite = false;
    } else {
        appendValue(mSeedTypes, newSeedType);
        appendValue(mMisMatchPos, pMisMatchPos);

        appendValue(mEffectiveSites, true);
        pEffectiveSite = true;
    }

    return pEffectiveSite;
}

} // namespace ts5cs
