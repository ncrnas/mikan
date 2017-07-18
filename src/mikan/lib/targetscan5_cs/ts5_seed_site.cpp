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

} // namespace ts5cs
