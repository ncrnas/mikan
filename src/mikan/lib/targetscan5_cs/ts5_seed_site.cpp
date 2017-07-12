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
    mBM = false;
    mLP = false;
    mOther = false;
    mAddInReverse = false;
}

//
// TS5SeedSites methods
//
void TS5SeedSites::reset_finder() {
    goBegin(mFinder);
    clear(mFinder);
}

int TS5SeedSites::find_seed_sites(mikan::TRNAStr const &pMiRNA) {
    TS5SeedSeqs seedSeqs;
    mikan::TRNAStr seedSeq;
    int retVal;

    reset_finder();

    seedSeqs.set_mirna_seq(pMiRNA);
    mikan::TCharSet mNullSet;
    seedSeqs.set_flags(mNullSet);
    retVal = seedSeqs.create_seed_seqs();
    if (retVal != 0) {
        std::cerr << "ERROR: Could not get the seed sequence for " << pMiRNA;
        std::cerr << "." << std::endl;
        return 1;
    }

    seedSeq = seedSeqs.get_seed_seq(0);
    while (find(mFinder, seedSeq)) {
        appendValue(mMRNAPos, position(mFinder).i1);
        appendValue(mSitePos, position(mFinder).i2);
    }

    return 0;
}

void TS5SeedSites::clear_pos() {
    clear(mMRNAPos);
    clear(mSitePos);
}

} // namespace ts5cs
