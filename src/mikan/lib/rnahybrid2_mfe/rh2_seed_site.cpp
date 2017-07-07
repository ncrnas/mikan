#include "mk_typedef.hpp"        // TRNATYPE
#include "rh2_seed_site.hpp"     // RH2SeedSeqs, RH2SeedSites
#include <seqan/seq_io.h>

using namespace seqan;
using namespace mikan;

namespace rh2mfe {

//
// RH2SeedSeqs methods
//
template<class TRNAString>
void RH2SeedSeqs<TRNAString>::set_mirna_seq(TRNAString pSeq) {
    clear(mSeedSeqs);
    clear(mSeedTypes);
    clear(mEffectiveSeeds);
    mMiRNASeq = pSeq;
}

template<class TRNAString>
int RH2SeedSeqs<TRNAString>::create_seed_seqs(CharString &pSeedDef, CharString &pOverlapDef) {
    if (length(mMiRNASeq) == 0) {
        return 1;
    }

    TRNAString seedSeq;
    CharString seedType;

    int retVal;

    resize(seedSeq, SEED_LEN);
    for (unsigned i = 0; i < length(seedSeq); ++i) {
        seedSeq[i] = mMiRNASeq[i + 1];
    }
    reverseComplement(seedSeq);

    if (pSeedDef[0] == '6' || pSeedDef[0] == '7') {
        (void) create_nmer_seed_seqs(seedSeq, pSeedDef);
    }

    resize(mEffectiveSeeds, length(mSeedSeqs), true);
    retVal = check_redundant_seeds(pOverlapDef);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

template<class TRNAString>
int RH2SeedSeqs<TRNAString>::create_nmer_seed_seqs(TRNAString &pSeedSeq, CharString &pSeedDef) {
    int retVal = 0;

    appendValue(mSeedSeqs, pSeedSeq);
    if (pSeedDef[0] == '6') {
        appendValue(mSeedTypes, "6mer");

    } else if (pSeedDef[0] == '7') {
        appendValue(mSeedTypes, "7mer");
    }

    if (pSeedDef[2] == 'G' && pSeedDef[3] == 'U') {
        CharString GUT;
        CharString GUM;

        if (pSeedDef[0] == '7') {
            GUT = "7mer_GUT";
            GUM = "7mer_GUM";
        } else if (pSeedDef[0] == '6') {
            GUT = "6mer_GUT";
            GUM = "6mer_GUM";
        } else {
            GUT = "GUT";
            GUM = "GUM";
        }

        if (pSeedDef == "6mGU+" || pSeedDef == "7mGU+") {
            retVal = create_multi_guwobble_seed_seqs(pSeedSeq, GUT, GUM);
        } else if (pSeedDef == "6mGU1" || pSeedDef == "7mGU1") {
            retVal = create_single_guwobble_seed_seqs(pSeedSeq, GUT, GUM);
        }
    }

    if (retVal != 0) {
        return 1;
    }

    return 0;
}

template<class TRNAString>
int RH2SeedSeqs<TRNAString>::create_single_guwobble_seed_seqs(
        TRNAString &pSeedSeq,
        seqan::CharString &pGUT,
        seqan::CharString &pGUM) {
    TRNAString seedGUSeq;

    for (unsigned i = 0; i < length(pSeedSeq); ++i) {
        if (pSeedSeq[i] == 'C') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'U';
            appendValue(mSeedSeqs, seedGUSeq);
            appendValue(mSeedTypes, pGUT);
        } else if (pSeedSeq[i] == 'A') {
            seedGUSeq = pSeedSeq;
            seedGUSeq[i] = 'G';
            appendValue(mSeedSeqs, seedGUSeq);
            appendValue(mSeedTypes, pGUM);
        }
    }

    return 0;
}

template<class TRNAString>
int RH2SeedSeqs<TRNAString>::create_multi_guwobble_seed_seqs(
        TRNAString &pSeedSeq,
        seqan::CharString &pGUT,
        seqan::CharString &pGUM) {
    TRNAString seedGUSeq;

    unsigned seedDatLen = (unsigned) length(mSeedSeqs);

    for (unsigned i = 0; i < length(pSeedSeq); ++i) {
        if (pSeedSeq[i] == 'C') {
            for (unsigned j = 0; j < seedDatLen; ++j) {
                seedGUSeq = get_seed_seq(j);
                seedGUSeq[i] = 'U';
                appendValue(mSeedSeqs, seedGUSeq);
                appendValue(mSeedTypes, pGUT);
            }
            seedDatLen = (unsigned) length(mSeedSeqs);
        } else if (pSeedSeq[i] == 'A') {
            for (unsigned j = 0; j < seedDatLen; ++j) {
                seedGUSeq = get_seed_seq(j);
                seedGUSeq[i] = 'G';
                appendValue(mSeedSeqs, seedGUSeq);
                appendValue(mSeedTypes, pGUM);
            }
            seedDatLen = (unsigned) length(mSeedSeqs);
        }

    }

    return 0;
}

template<class TRNAString>
int RH2SeedSeqs<TRNAString>::check_redundant_seeds(seqan::CharString &) {
    typedef Index<StringSet<TRNAString>, IndexQGram<UngappedShape<SEED_LEN> > > TIndexQGram;
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
// RH2SeedSites methods
//
template<class TRNAString>
void RH2SeedSites<TRNAString>::reset_finder() {
    goBegin(mFinder);
    clear(mFinder);
}

template<class TRNAString>
int RH2SeedSites<TRNAString>::find_seed_sites(
        TRNAString const &pMiRNA,
        CharString &pSeedDef,
        CharString &pOverlapDef) {
    RH2SeedSeqs<TRNAString> seedSeqs;
    TRNAString seedSeq;
    CharString seedType;
    int retVal;
    unsigned mRNAPos, sitePos;
    bool effectiveSite;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    retVal = seedSeqs.create_seed_seqs(pSeedDef, pOverlapDef);
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
            set_new_seed_type(seedType, pSeedDef, mRNAPos, sitePos, pMiRNA, effectiveSite);
            appendValue(mEffectiveSites, effectiveSite);
        }
        reset_finder();
    }

    return 0;
}

template<class TRNAString>
void RH2SeedSites<TRNAString>::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mEffectiveSites);
}

template<class TRNAString>
void RH2SeedSites<TRNAString>::set_new_seed_type(
        CharString &pCurSeedType,
        CharString &pSeedDef,
        unsigned pMRNAPos,
        unsigned pSitePos,
        TRNAString const &pMiRNA,
        bool &pEffectiveSite) {
    bool IsA1, MatchM8;
    TRNAString revMiRNASeq, miRNAM8, mRNAM8, miRNAM8c, mRNAA1;
    CharString newSeedType = "";
    TRNAString miRNASeq = pMiRNA;
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

// Explicit template instantiation

template
class RH2SeedSeqs<TRNATYPE>;

template
class RH2SeedSites<TRNATYPE>;

} // namespace rh2mfe
