#include <seqan/seq_io.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "ts5_seed_site.hpp"     // TS5SeedSeqs, TS5SeedSites

using namespace seqan;

namespace ts5cs {

//
// TS5SeedSeqs methods
//
int TS5SeedSeqs::create_seed_seqs() {
    if (length(mMiRNASeq) == 0) {
        return 1;
    }

    mikan::TRNAStr seedSeq;

    resize(seedSeq, 6);
    for (unsigned i = 0; i < length(seedSeq); ++i) {
        seedSeq[i] = mMiRNASeq[i + 1];
    }
    reverseComplement(seedSeq);

    appendValue(mSeedSeqs, seedSeq);

    return 0;
}

void TS5SeedSeqs::set_mirna_seq(mikan::TRNAStr pSeq) {
    clear(mSeedSeqs);
    mMiRNASeq = pSeq;
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
    retVal = seedSeqs.create_seed_seqs();
    if (retVal != 0) {
        std::cerr << "ERROR: Could not get the seed sequence for " << pMiRNA;
        std::cerr << "." << std::endl;
        return 1;
    }

    seedSeq = seedSeqs.get_seed_seq();
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
