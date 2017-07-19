#include <seqan/seq_io.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_seed_site.hpp"      // TSSVMSequences, TSSVMSeedSeqs, TSSVMSeedSites, TSSVMSeedSiteOverlap

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
// TSSVMSeedSiteOverlap methods
//
void TSSVMSeedSiteOverlap::clear_site_pos() {
    mRNAPosSet.clear();
    mSiteMap.clear();
    clear(mSortedMRNAPos);
}

void TSSVMSeedSiteOverlap::cluster_site_pos(TSSVMSeedSites &pSeedSites) {
    const String<unsigned> &mRNAPos = pSeedSites.get_mrna_pos();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        mRNAPosSet.insert(mRNAPos[i]);
        mSiteMap.insert(std::pair<unsigned, unsigned>(mRNAPos[i], i));
    }
}

int TSSVMSeedSiteOverlap::filter_overlapped_sites(
        TSSVMSeedSites &pSeedSites,
        unsigned pMRNANum) {
    TItSet itSet;

    cluster_site_pos(pSeedSites);
    resize(mSortedMRNAPos, pMRNANum);

    for (itSet = mRNAPosSet.begin(); itSet != mRNAPosSet.end(); ++itSet) {
        sort_by_seed_type(pSeedSites, *itSet);
    }

    return 0;
}

void TSSVMSeedSiteOverlap::sort_by_seed_type(
        TSSVMSeedSites &pSeedSites,
        int pPosIdx) {
    TItMap itMap;
    TItRetPair ret;
    std::multimap<unsigned, unsigned> sortedPos;
    std::multimap<unsigned, unsigned> sortedSeeds;
    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    unsigned preced;
    unsigned count = 0;
    TITPos itPos;

    ret = mSiteMap.equal_range(pPosIdx);
    for (itMap = ret.first; itMap != ret.second; ++itMap) {
        sortedPos.insert(TPosPair((unsigned) sitePos[(*itMap).second], (*itMap).second));
        ++count;
    }

    resize(mSortedMRNAPos[pPosIdx], count);
    int n = 0;
    for (itPos = sortedPos.begin(); itPos != sortedPos.end(); ++itPos) {
        preced = get_seedtype_precedence(seedTypes[(*itPos).second]);
        sortedSeeds.insert(TPosPair(preced * (count + 1) + n, (*itPos).second));
        mSortedMRNAPos[pPosIdx][n] = (*itPos).second;
        ++n;
    }

    if (count > 1) {
        mark_overlapped_sites(pSeedSites, sortedSeeds);
    }
}

unsigned TSSVMSeedSiteOverlap::get_seedtype_precedence(const CharString &pSeedType) {
    unsigned preced;

    if (pSeedType == "8mer" || pSeedType == "7mer-A1" || pSeedType == "7mer-m8") {
        preced = 0;
    } else if (pSeedType == "6mer") {
        preced = 1;
    } else if (pSeedType == "GUM") {
        preced = 2;
    } else if (pSeedType == "GUT") {
        preced = 3;
    } else if (pSeedType == "BT") {
        preced = 4;
    } else if (pSeedType == "BM") {
        preced = 5;
    } else if (pSeedType == "LP") {
        preced = 6;
    } else {
        preced = 7;
    }

    return preced;
}

void TSSVMSeedSiteOverlap::mark_overlapped_sites(
        TSSVMSeedSites &pSeedSites,
        std::multimap<unsigned, unsigned> &pSortedSeeds) {
    IntervalTree<unsigned> tree;
    String<unsigned> results;
    const String<unsigned> &sitePos = pSeedSites.get_site_pos();
    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    TITSeedTypes itSeedType;
    unsigned pos7, pos1, startPos, endPos;
    CharString seedType;
    bool searchOverlap;

    startPos = 0;
    endPos = 0;
    for (itSeedType = pSortedSeeds.begin(); itSeedType != pSortedSeeds.end(); ++itSeedType) {
        seedType = seedTypes[(*itSeedType).second];

        pos7 = sitePos[(*itSeedType).second] - 1;
        pos1 = sitePos[(*itSeedType).second] + 7;
        searchOverlap = true;

        if (seedType == "8mer" || seedType == "7mer-A1" || seedType == "7mer-m8" || seedType == "6mer") {
            startPos = pos7;
            endPos = pos1;
            searchOverlap = false;
        } else if (seedType == "GUM" || seedType == "GUT" || seedType == "LP") {
            startPos = pos7 - 1;
            endPos = pos1 + 1;
        } else if (seedType == "BT") {
            startPos = pos7 - 1;
            endPos = pos1 + 2;
        } else if (seedType == "BM") {
            pos7 = sitePos[(*itSeedType).second];
            pos1 = sitePos[(*itSeedType).second] + 8;
            startPos = pos7 - 1;
            endPos = pos1;
        }

        clear(results);
        if (searchOverlap) {
            findIntervals(tree, startPos, endPos, results);
        }

        if (length(results) == 0) {
            addInterval(tree, pos7, pos1);
        }

        if (searchOverlap && length(results) > 0) {
            pSeedSites.mEffectiveSites[(*itSeedType).second] = false;
        }
    }
}

} // namespace tssvm
