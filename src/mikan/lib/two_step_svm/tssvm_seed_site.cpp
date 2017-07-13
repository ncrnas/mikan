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
void TSSVMSeedSites::reset_finder() {
    goBegin(mFinder);
    clear(mFinder);
}

int TSSVMSeedSites::find_seed_sites(
        mikan::TRNAStr const &pMiRNA) {
    TSSVMSeedSeqs seedSeqs;
    mikan::TRNAStr seedSeq;
    CharString newSeedType;
    int retVal;
    unsigned mRNAPos, sitePos;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    mikan::TCharSet mNullSet;
    seedSeqs.set_flags(mNullSet);
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

        while (find(mFinder, seedSeq)) {
            mRNAPos = position(mFinder).i1;
            sitePos = position(mFinder).i2;

            const CharString &seedType = seedSeqs.get_seed_type(i);
            if (seedType == "GUM" || seedType == "GUT" || seedType == "LP" || seedType == "BT") {
                if ((sitePos < MIN_DIST_TO_CDS + 1) ||
                    (sitePos + MIN_DIST_UTR_END > length(mMRNASeqs[mRNAPos]) - 1)) {
                    continue;
                }
            } else {
                if ((sitePos < MIN_DIST_TO_CDS) || (sitePos + MIN_DIST_UTR_END > length(mMRNASeqs[mRNAPos]))) {
                    continue;
                }
            }

            retVal = set_seed_pos(seedSeqs, mRNAPos, sitePos, pMiRNA, seedSeq, i);
            if (retVal == 0) {
                appendValue(mMRNAPos, mRNAPos);
                appendValue(mSitePos, sitePos);
                appendValue(mEffectiveSites, true);
            }
        }
        reset_finder();
    }

    return 0;
}

void TSSVMSeedSites::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mS1Pos);
    clear(mS8Pos);
    clear(mMisMatchPos);
    clear(mEffectiveSites);
}

int TSSVMSeedSites::set_seed_pos(
        TSSVMSeedSeqs &pSeedSeqs,
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr const &pMiRNA,
        const mikan::TRNAStr &pSeedSeq,
        unsigned pIdx) {
    unsigned m8Pos, a1Pos;
    const CharString &seedType = pSeedSeqs.get_seed_type(pIdx);
    unsigned mmPos = pSeedSeqs.get_mismatched_pos(pIdx);
    int retVal;

    m8Pos = pSitePos - 1;
    if (seedType == "BM") {
        m8Pos = pSitePos;
    }

    a1Pos = pSitePos + 6;
    if (seedType == "BT") {
        a1Pos = pSitePos + 7;
    }

    if ((seedType != "6mer") && (a1Pos >= length(mMRNASeqs[pMRNAPos]))) {
        return 1;
    }

    retVal = set_seed_type(seedType, mMRNASeqs[pMRNAPos], pMiRNA, m8Pos, a1Pos, mmPos, pSeedSeq);
    if (retVal != 0) {
        return 1;
    }

    appendValue(mS1Pos, a1Pos);
    appendValue(mS8Pos, m8Pos);

    return 0;

}

int TSSVMSeedSites::set_seed_type(
        const CharString &pCurType,
        const mikan::TRNAStr &pMRNASeq,
        const mikan::TRNAStr &pMiRNASeq,
        unsigned pM8Pos,
        unsigned pA1Pos,
        unsigned pMisMatchedPos,
        const mikan::TRNAStr &) {
    mikan::TRNAStr newSeedSeq, complMiRNASeq;
    seqan::Rna miRNAM2, miRNAM8;
    seqan::Rna mRNAM2, mRNAM8, mRNAA1;
    CharString newSeedType;
    bool IsA1, MatchM8;

    complMiRNASeq = pMiRNASeq;
    complement(complMiRNASeq);
    miRNAM8 = complMiRNASeq[7];
    miRNAM2 = complMiRNASeq[1];

    mRNAM8 = pMRNASeq[pM8Pos];
    mRNAM2 = pMRNASeq[pA1Pos - 1];
    if (pA1Pos < length(pMRNASeq)) {
        mRNAA1 = pMRNASeq[pA1Pos];
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
    if (pCurType == "6mer") {
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
        if (pCurType == "BT") {
            if (miRNAM2 == mRNAM2) {
                newSeedType = pCurType;
            }
        } else {
            newSeedType = pCurType;
        }
    }

    if (newSeedType == "") {
        return 1;
    }

    appendValue(mSeedTypes, newSeedType);
    appendValue(mMisMatchPos, pMisMatchedPos);

    return 0;
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
