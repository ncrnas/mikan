#include <seqan/seq_io.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"       // MKSeedSeqs
#include "mk_seed_site.hpp"      // MKSeedSites

using namespace seqan;

namespace mikan {

//
// MKSeedSites methods
//
void MKSeedSites::reset_finder() {
    goBegin(mFinder);
    clear(mFinder);
}

void MKSeedSites::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
    clear(mSeedTypes);
    clear(mMisMatchPos);
    clear(mEffectiveSites);

    clear(mS1Pos);
    clear(mS8Pos);

}

int MKSeedSites::find_seed_sites(mikan::MKSeedSeqs &seedSeqs, mikan::TCharSet &pSeedTypeDef) {
    mikan::TRNAStr seedSeq;
    CharString seedType;
    unsigned mRNAPos, sitePos;
    bool effectiveSite, effectiveSite2;
    unsigned misMatchPos;

    mikan::TRNAStr miRNASeq = seedSeqs.get_mirna_seq();

    reset_finder();

    for (unsigned i = 0; i < length(seedSeqs.mEffectiveSeeds); ++i) {
        if (!seedSeqs.mEffectiveSeeds[i]) {
            continue;
        }

        seedSeq = seedSeqs.get_seed_seq(i);
        seedType = seedSeqs.get_seed_type(i);
        misMatchPos = seedSeqs.get_mismatched_pos(i);

        while (find(mFinder, seedSeq)) {
            mRNAPos = position(mFinder).i1;
            sitePos = position(mFinder).i2;

            effectiveSite = check_position_1(mRNAPos, sitePos, seedType);
            if (effectiveSite) {
                effectiveSite = set_new_seed_type(mRNAPos, sitePos, miRNASeq, pSeedTypeDef, seedType, misMatchPos,
                                                  effectiveSite);
            }
            if (effectiveSite) {
                appendValue(mMRNAPos, mRNAPos);
                appendValue(mSitePos, sitePos);
                effectiveSite2 = check_position_2(mRNAPos, sitePos, seedType);
                if (!effectiveSite2) {
                    unsigned curIdx = length(mEffectiveSites) - 1;
                    mEffectiveSites[curIdx] = false;

                }
            }
        }
        reset_finder();
    }

    if (mUpdatePos) {
        resize(mS1Pos, length(mEffectiveSites));
        resize(mS8Pos, length(mEffectiveSites));

        for (unsigned i = 0; i < length(mEffectiveSites); i++) {
            mS1Pos[i] = mSitePos[i] + 6;
            mS8Pos[i] = mSitePos[i] - 1;
        }

    }

    return 0;
}

bool MKSeedSites::check_position_1(unsigned pMRNAPos, unsigned pSitePos, seqan::CharString &) {
    bool effectiveSite = true;

    if ((pSitePos < mMinToCDS) || (pSitePos + mMinToEnd > length(mMRNASeqs[pMRNAPos]))) {
        effectiveSite = false;;
    }

    return effectiveSite;
}

bool MKSeedSites::set_new_seed_type(
        unsigned,
        unsigned,
        mikan::TRNAStr &,
        mikan::TCharSet &,
        seqan::CharString &pSeedType,
        int,
        bool pEffectiveSite) {

    appendValue(mSeedTypes, pSeedType);
    appendValue(mMisMatchPos, 0);
    appendValue(mEffectiveSites, pEffectiveSite);

    return pEffectiveSite;
}

void MKSeedSites::set_mx_matches(
        unsigned pMRNAPos,
        unsigned pSitePos,
        mikan::TRNAStr const &pMiRNA,
        int pMx,
        bool &pNoMx,
        bool &pMatchMx,
        bool &pGutMx,
        bool &pGumMx,
        bool &pIsA) {

    int pos;
    mikan::TRNAStr cMiRNASeq, miRNAMx, mRNAMx, miRNAMxC;

    miRNAMx = pMiRNA[pMx - 1];
    cMiRNASeq = pMiRNA;
    complement(cMiRNASeq);
    miRNAMxC = cMiRNASeq[pMx - 1];

    pMatchMx = false;
    pGutMx = false;
    pGumMx = false;
    pIsA = false;

    pos = static_cast<int>(pSitePos) - pMx + 1 + static_cast<int>(INDEXED_SEQ_LEN);
    if (pos >= 0 && pos < static_cast<int>(length(mMRNASeqs[pMRNAPos]))) {
        mRNAMx = mMRNASeqs[pMRNAPos][pos];
        if (mRNAMx == 'A') {
            pIsA = true;
        }
        if (miRNAMxC == mRNAMx) {
            pMatchMx = true;
        } else if (miRNAMx == 'G' && mRNAMx == 'U') {
            pGutMx = true;
        } else if (miRNAMx == 'U' && mRNAMx == 'G') {
            pGumMx = true;
        }
        pNoMx = false;
    } else {
        pNoMx = true;
    }
}

void MKSeedSites::print_all() {
    for (unsigned i = 0; i < length(mEffectiveSites); i++) {
        std::cout << i << ", ";
        std::cout << mSeedTypes[i] << ", ";
        std::cout << mMRNAPos[i] << ", ";
        std::cout << mSitePos[i] << ", ";
        std::cout << mMisMatchPos[i] << ", ";
        std::cout << mEffectiveSites[i] << std::endl;
    }
}

} // namespace mikan
