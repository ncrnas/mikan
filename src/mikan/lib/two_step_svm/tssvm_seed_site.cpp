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
bool TSSVMSeedSites::check_position(unsigned pMRNAPos, unsigned pSitePos, seqan::CharString &pSeedType) {
    bool effectiveSite;

    effectiveSite = true;
    if (pSeedType == "GUM" || pSeedType == "GUT" || pSeedType == "LP" || pSeedType == "BT") {
        if ((pSitePos < MIN_DIST_TO_CDS + 1) ||
            (pSitePos + MIN_DIST_UTR_END > length(mMRNASeqs[pMRNAPos]) - 1)) {
            effectiveSite = false;;
        }
    } else {
        if ((pSitePos < MIN_DIST_TO_CDS) || (pSitePos + MIN_DIST_UTR_END > length(mMRNASeqs[pMRNAPos]))) {
            effectiveSite = false;;
        }
    }

    return effectiveSite;
}

void TSSVMSeedSites::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mMisMatchPos);
    clear(mEffectiveSites);

    clear(mS1Pos);
    clear(mS8Pos);
}

bool TSSVMSeedSites::set_new_seed_type(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr &pMiRNASeq,
        mikan::TCharSet &,
        seqan::CharString &pSeedType,
        int pMisMatchPos,
        bool pEffectiveSite) {

    unsigned m8Pos, a1Pos;
    CharString newSeedType;

    mikan::TRNAStr newSeedSeq, complMiRNASeq;
    seqan::Rna miRNAM2, miRNAM8;
    seqan::Rna mRNAM2, mRNAM8, mRNAA1;
    bool IsA1, MatchM8;

    m8Pos = pSitePos - 1;
    if (pSeedType == "BM") {
        m8Pos = pSitePos;
    }

    a1Pos = pSitePos + 6;
    if (pSeedType == "BT") {
        a1Pos = pSitePos + 7;
    }
    if ((pSeedType != "6mer") && (a1Pos >= length(mMRNASeqs[pMRNAPos]))) {
        return false;
    }

    complMiRNASeq = pMiRNASeq;
    complement(complMiRNASeq);
    miRNAM8 = complMiRNASeq[7];
    miRNAM2 = complMiRNASeq[1];

    mRNAM8 = mMRNASeqs[pMRNAPos][m8Pos];
    mRNAM2 = mMRNASeqs[pMRNAPos][a1Pos - 1];
    if (a1Pos < length(mMRNASeqs[pMRNAPos])) {
        mRNAA1 = mMRNASeqs[pMRNAPos][a1Pos];
    } else {
        mRNAA1 = 'A';
    }

    MatchM8 = false;
    if (miRNAM8 == mRNAM8) {
        MatchM8 = true;
    }

    IsA1 = false;
    if (mRNAA1 == 'A') {
        IsA1 = true;
    }

    newSeedType = "";
    if (pSeedType == "6mer") {
        if (IsA1 && MatchM8) {
            newSeedType = "8mer";
        } else if (IsA1) {
            newSeedType = "7mer-A1";
        } else if (MatchM8) {
            newSeedType = "7mer-m8";
        } else {
            newSeedType = "6mer";
        }
    } else if (IsA1 && MatchM8) {
        if (pSeedType == "BT") {
            if (miRNAM2 == mRNAM2) {
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

        appendValue(mS1Pos, a1Pos);
        appendValue(mS8Pos, m8Pos);
        
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
