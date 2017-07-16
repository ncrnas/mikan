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
}

int MKSeedSites::find_seed_sites(mikan::MKSeedSeqs &seedSeqs, mikan::TCharSet &pSeedTypeDef) {
    mikan::TRNAStr seedSeq;
    CharString seedType;
    unsigned mRNAPos, sitePos;
    bool effectiveSite;
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

            effectiveSite = check_position(mRNAPos, sitePos);

            appendValue(mMRNAPos, mRNAPos);
            appendValue(mSitePos, sitePos);

            if (effectiveSite) {
                set_new_seed_type(mRNAPos, sitePos, miRNASeq, pSeedTypeDef, seedType, misMatchPos, effectiveSite);
            } else {
                appendValue(mSeedTypes, "");
                appendValue(mMisMatchPos, 0);
                appendValue(mEffectiveSites, false);
            }
        }
        reset_finder();
    }

    return 0;
}

bool MKSeedSites::check_position(unsigned, unsigned) {
    return true;
}

void MKSeedSites::set_new_seed_type(
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

}

void MKSeedSites::print_all() {
    for (unsigned i = 0 ; i < length(mEffectiveSites); i++) {
        std::cout << i << ", ";
        std::cout << mSeedTypes[i] << ", ";
        std::cout << mMRNAPos[i] << ", ";
        std::cout << mSitePos[i] << ", ";
        std::cout << mMisMatchPos[i] << ", ";
        std::cout << mEffectiveSites[i] << std::endl;
    }
}

} // namespace mikan
