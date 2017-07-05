#include <tm1_inst_template.hpp> // TRNATYPE
#include <tm1_seed_site.hpp>     // TM1Sequences, TM1SeedSeqs, TM1SeedSites
#include <seqan/seq_io.h>

using namespace seqan;

namespace tm1p {

//
// TM1Sequences methods
//
template<class TRNAString>
int TM1Sequences<TRNAString>::read_fasta(CharString const &pFasta) {
    CharString id;
    CharString seq;

    SequenceStream seqStream(toCString(pFasta));
    if (!isGood(seqStream)) {
        std::cerr << "ERROR: Could not open the file!" << std::endl;
        return 1;
    }

    while (!atEnd(seqStream)) {
        if (readRecord(id, seq, seqStream) != 0) {
            std::cerr << "ERROR: Could not read from " << toCString(pFasta) << "!" << std::endl;
            return 1;
        }

        toUpper(seq);
        for (unsigned i = 0; i < length(seq); ++i) {
            if (seq[i] == 'T') {
                seq[i] = 'U';
            }
        }

        if (length(seq) > 0) {
            appendValue(mSeqIds, id);
            appendValue(mSeqs, seq);
        }
    }

    return 0;
}

//
// TM1SeedSeqs methods
//
template<class TRNAString>
void TM1SeedSeqs<TRNAString>::set_mirna_seq(TRNAString pSeq) {
    clear(mSeedSeqs);
    clear(mSeedTypes);
    clear(mEffectiveSeeds);
    mMiRNASeq = pSeq;
}

template<class TRNAString>
int TM1SeedSeqs<TRNAString>::create_seed_seqs() {
    if (length(mMiRNASeq) == 0) {
        return 1;
    }

    TRNAString seedSeq;

    int retVal;

    resize(seedSeq, 6);
    for (unsigned i = 0; i < length(seedSeq); ++i) {
        seedSeq[i] = mMiRNASeq[i + 1];
    }
    reverseComplement(seedSeq);

    (void) create_nmer_seed_seqs(seedSeq);

    resize(mEffectiveSeeds, length(mSeedSeqs), true);
    retVal = check_redundant_seeds();
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

template<class TRNAString>
int TM1SeedSeqs<TRNAString>::create_nmer_seed_seqs(TRNAString &pSeedSeq) {
    int retVal = 0;

    appendValue(mSeedSeqs, pSeedSeq);
    appendValue(mSeedTypes, "6mer");

    CharString lpSeedType = "MM";
    retVal = create_lp_seed_seqs(pSeedSeq, lpSeedType);
    if (retVal != 0) {
        return 1;
    }

    retVal = create_single_guwobble_seed_seqs(pSeedSeq);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

template<class TRNAString>
int TM1SeedSeqs<TRNAString>::create_single_guwobble_seed_seqs(TRNAString &pSeedSeq) {
    TRNAString seedGUSeq;
    int retVal;
    CharString lpSeedType = "GUMM";

    for (unsigned i = 0; i < length(pSeedSeq); ++i) {
        if (pSeedSeq[i] == 'C') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'U';
            appendValue(mSeedSeqs, seedGUSeq);
            appendValue(mSeedTypes, "GUT");
            retVal = create_lp_seed_seqs(seedGUSeq, lpSeedType);
            if (retVal != 0) {
                return 1;
            }
        } else if (pSeedSeq[i] == 'A') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'G';
            appendValue(mSeedSeqs, seedGUSeq);
            appendValue(mSeedTypes, "GUM");
            retVal = create_lp_seed_seqs(seedGUSeq, lpSeedType);
            if (retVal != 0) {
                return 1;
            }
        }
    }

    return 0;
}

template<class TRNAString>
int TM1SeedSeqs<TRNAString>::create_lp_seed_seqs(TRNAString &pSeedSeq, CharString &pSeedType) {
    TRNAString seedLPSeq;
    char ch1 = 0;
    char ch2 = 0;
    char ch3 = 0;
    int m2pos = length(pSeedSeq) - 1;
    CharString seedType2;

    seedType2 = pSeedType;
    if (pSeedSeq[m2pos] == 'A') {
        ch1 = 'C';
        ch2 = 'G';
        ch3 = 'U';
        if (pSeedType == "GUMM") {
            seedType2 = "GUGU";
        } else {
            seedType2 = "GUM";
        }
    } else if (pSeedSeq[m2pos] == 'C') {
        ch1 = 'A';
        ch2 = 'U';
        ch3 = 'G';
        if (pSeedType == "GUMM") {
            seedType2 = "GUGU";
        } else {
            seedType2 = "GUT";
        }
    } else if (pSeedSeq[m2pos] == 'G') {
        ch1 = 'A';
        ch2 = 'U';
        ch3 = 'C';
    } else if (pSeedSeq[m2pos] == 'U') {
        ch1 = 'C';
        ch2 = 'G';
        ch3 = 'A';
    }

    seedLPSeq = pSeedSeq;
    seedLPSeq[m2pos] = ch1;
    appendValue(mSeedSeqs, seedLPSeq);
    appendValue(mSeedTypes, pSeedType);

    seedLPSeq = pSeedSeq;
    seedLPSeq[m2pos] = ch2;
    appendValue(mSeedSeqs, seedLPSeq);
    appendValue(mSeedTypes, seedType2);

    seedLPSeq = pSeedSeq;
    seedLPSeq[m2pos] = ch3;
    appendValue(mSeedSeqs, seedLPSeq);
    appendValue(mSeedTypes, pSeedType);

    return 0;
}

template<class TRNAString>
int TM1SeedSeqs<TRNAString>::check_redundant_seeds() {
    typedef Index<StringSet<TRNAString>, IndexQGram<UngappedShape<6> > > TIndexQGram;
    typedef Finder<TIndexQGram> TFinder;

    TRNAString seedSeq;
    TIndexQGram RNAIdx(mSeedSeqs);
    TFinder finder(RNAIdx);

    for (unsigned i = 0; i < length(mSeedSeqs); ++i) {
        if (!mEffectiveSeeds[i]) {
            continue;
        }
        seedSeq = mSeedSeqs[i];
        while (find(finder, seedSeq)) {
            if (i != position(finder).i1) {
                mEffectiveSeeds[position(finder).i1] = false;
            }
        }
        goBegin(finder);
        clear(finder);
    }

    return 0;
}

//
// TM1SeedSites methods
//
template<class TRNAString>
void TM1SeedSites<TRNAString>::reset_finder() {
    goBegin(mFinder);
    clear(mFinder);
}

template<class TRNAString>
int TM1SeedSites<TRNAString>::find_seed_sites(TRNAString const &pMiRNA) {
    TM1SeedSeqs<TRNAString> seedSeqs;
    TRNAString seedSeq;
    CharString seedType;
    int retVal;
    unsigned mRNAPos, sitePos;
    bool effectiveSite;
    unsigned idx;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    retVal = seedSeqs.create_seed_seqs();
    if (retVal != 0) {
        std::cerr << "ERROR: Could not get the seed sequence for " << pMiRNA;
        std::cerr << std::endl;
        return 1;
    }

    idx = 0;
    for (unsigned i = 0; i < length(seedSeqs.mEffectiveSeeds); ++i) {
        if (!seedSeqs.mEffectiveSeeds[i]) {
            continue;
        }

        seedSeq = seedSeqs.get_seed_seq(i);
        seedType = seedSeqs.get_seed_type(i);

        while (find(mFinder, seedSeq)) {
            mRNAPos = position(mFinder).i1;
            sitePos = position(mFinder).i2;

            appendValue(mMRNAPos, mRNAPos);
            appendValue(mSitePos, sitePos);
            effectiveSite = false;
            set_new_seed_type(seedType, mRNAPos, sitePos, pMiRNA, effectiveSite);
            appendValue(mEffectiveSites, effectiveSite);
            ++idx;
        }
        reset_finder();
    }

    return 0;
}

template<class TRNAString>
void TM1SeedSites<TRNAString>::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mEffectiveSites);
    clear(mM8Match);
    clear(mM8GU);
    clear(mM1A);
    clear(mM1Match);
    clear(mM1GU);
    clear(mMRNASeqLen);
    clear(mM8Pos);
}

template<class TRNAString>
void TM1SeedSites<TRNAString>::get_mx_match(
        TRNAString const &pMiRNASeq,
        TRNAString const &pMiRNACompSeq,
        TRNAString const &pMRNASeq,
        unsigned pSitePos,
        int pMx,
        bool &pMatch,
        bool &pGU,
        bool &pNoMx,
        bool &pIsA) {
    int pos;
    TRNAString miRNAMx, miRNAMxC, mRNAMx;

    pos = pSitePos - (pMx - 7);
    miRNAMx = pMiRNASeq[pMx - 1];
    miRNAMxC = pMiRNACompSeq[pMx - 1];

    pMatch = false;
    pGU = false;
    pIsA = false;

    if (pos >= 0 && pos < (int) length(pMRNASeq)) {
        mRNAMx = pMRNASeq[pos];
        if (mRNAMx == 'A') {
            pIsA = true;
        }
        if (miRNAMxC == mRNAMx) {
            pMatch = true;
        }
        if (((miRNAMx == 'G') && (mRNAMx == 'U')) || ((miRNAMx == 'U') && (mRNAMx == 'G'))) {
            pGU = true;
        }
        pNoMx = false;
    } else {
        pNoMx = true;
    }
}

template<class TRNAString>
void TM1SeedSites<TRNAString>::get_match_count(
        TRNAString const &pMiRNASeq,
        TRNAString const &pMiRNACompSeq,
        TRNAString const &pMRNASeq,
        unsigned pSitePos,
        int pMx1,
        int pMx2,
        int &pMatchCount,
        int &pGUCount) {
    bool isAx, matchMx, guMx, noMx;
    pMatchCount = 0;
    pGUCount = 0;

    for (int i = pMx1; i < pMx2 + 1; ++i) {
        get_mx_match(pMiRNASeq, pMiRNACompSeq, pMRNASeq, pSitePos, i, matchMx, guMx, noMx, isAx);
        if (!noMx) {
            if (matchMx) {
                ++pMatchCount;
            } else if (guMx) {
                ++pGUCount;
            }
        }
    }
}

template<class TRNAString>
void TM1SeedSites<TRNAString>::set_new_seed_type(
        CharString &,
        unsigned pMRNAPos,
        unsigned pSitePos,
        TRNAString const &pMiRNA,
        bool &pEffectiveSite) {
    CharString newSeedType = "";
    bool isA1, isAx;
    bool matchM1, matchM2, matchM8, matchM9;
    bool guM1, guM2, guM8, guM9;
    bool noM1, noM2, noM8, noM9;
    int matchCount, guCount, matchCount2, guCount2;
    TRNAString revMiRNASeq;

    revMiRNASeq = pMiRNA;
    complement(revMiRNASeq);

    // check m8, m2, m1 match
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 9, matchM9, guM9, noM9, isAx);
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 8, matchM8, guM8, noM8, isAx);
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 2, matchM2, guM2, noM2, isAx);
    get_mx_match(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 1, matchM1, guM1, noM1, isA1);

    get_match_count(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 2, 7, matchCount, guCount);
    get_match_count(pMiRNA, revMiRNASeq, mMRNASeqs[pMRNAPos], pSitePos, 2, 8, matchCount2, guCount2);

    if (isA1 && matchM8 && (matchCount + guCount == 6) && guCount < 2) {
        newSeedType = "8mer";
    } else if (!isA1 && matchM8 && (matchCount + guCount == 6) && guCount < 2) {
        newSeedType = "7mer-m8";
    } else if (isA1 && !matchM8 && (matchCount + guCount == 6) && guCount < 2) {
        newSeedType = "7mer-A1";
    } else if (((matchCount2 + guCount2 == 6) && guCount2 < 2) || ((matchCount2 + guCount2 > 6) && guCount2 < 3)) {
        newSeedType = "6mer";
    }

    if (newSeedType != "") {
        pEffectiveSite = true;
    }

    appendValue(mSeedTypes, newSeedType);

    appendValue(mM8Match, matchM8);
    appendValue(mM8GU, guM8);
    appendValue(mM1A, isA1);
    appendValue(mM1Match, matchM1);
    appendValue(mM1GU, guM1);

    appendValue(mMRNASeqLen, length(mMRNASeqs[pMRNAPos]));
    appendValue(mM8Pos, pSitePos - 1);
}

template<class TRNAString>
int TM1SeedSites<TRNAString>::get_seed_len(int pIdx) {
    if (mSeedTypes[pIdx] == "8mer") {
        return 8;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        return 7;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        return 7;
    } else if (mSeedTypes[pIdx] == "6mer") {
        return 6;
    }

    return 6;
}

template<class TRNAString>
int TM1SeedSites<TRNAString>::get_seed_start_pos(int pIdx) {
    int offset = 0;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 0;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 0;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match_gu(pIdx)) {
            offset = 0;
        } else {
            offset = 1;
        }
    }

    return std::max((int) mM8Pos[pIdx] + offset, 0);
}

template<class TRNAString>
int TM1SeedSites<TRNAString>::get_seed_start_pos2(int pIdx) {
    int offset = 1;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match_gu(pIdx)) {
            offset = 0;
        } else {
            offset = 1;
        }
    }

    return std::max((int) mM8Pos[pIdx] + offset, 0);
}

template<class TRNAString>
int TM1SeedSites<TRNAString>::get_length_to_cds(int pIdx) {
    int offset = 0;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 0;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 1;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match(pIdx)) {
            offset = 0;
        } else {
            offset = 1;
        }
    }

    return std::max((int) mM8Pos[pIdx] + offset, 0);
}

template<class TRNAString>
int TM1SeedSites<TRNAString>::get_seed_end_pos(int pIdx) {
    int offset = 0;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 7;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match_gu(pIdx)) {
            offset = 6;
        } else {
            offset = 7;
        }
    }


    return std::min((int) mM8Pos[pIdx] + offset, (int) mMRNASeqLen[pIdx]);
}


template<class TRNAString>
int TM1SeedSites<TRNAString>::get_seed_end_pos2(int pIdx) {
    int offset = 0;

    if (mSeedTypes[pIdx] == "8mer") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-A1") {
        offset = 8;
    } else if (mSeedTypes[pIdx] == "7mer-m8") {
        offset = 7;
    } else if (mSeedTypes[pIdx] == "6mer") {
        if (is_m8_match(pIdx)) {
            offset = 6;
        } else {
            offset = 7;
        }
    }


    return std::min((int) mM8Pos[pIdx] + offset, (int) mMRNASeqLen[pIdx]);
}

// Explicit template instantiation
template
class TM1Sequences<TRNATYPE>;

template
class TM1SeedSeqs<TRNATYPE>;

template
class TM1SeedSites<TRNATYPE>;

} // namespace tm1p
