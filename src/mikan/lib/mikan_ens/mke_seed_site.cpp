#include <seqan/seq_io.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mke_seed_site.hpp"     // MKESeedSeqs, MKESeedSites

using namespace seqan;

namespace mkens {

//
// MKESeedSeqs methods
//
void MKESeedSeqs::set_flags() {
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

} // namespace mkens

