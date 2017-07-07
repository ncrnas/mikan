#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_input.hpp"          // MKInput

namespace mikan {

//
// MKInput methods
//
void MKInput::set_file_names(seqan::CharString &pMiRNAFasta,
                                         seqan::CharString &pMRNAFasta) {
    mMiRNAFasta = pMiRNAFasta;
    mMRNAFasta = pMRNAFasta;
}

void MKInput::set_options(MKOptions & opt) {
    mMiRNAFasta = opt.mMiRNAFasta;
    mMRNAFasta = opt.mMRNAFasta;
}

int MKInput::load_seq_from_file() {
    int retVal;

    // Read miRNA fasta file
    retVal = mMiRNASeqs.read_fasta(mMiRNAFasta);
    if (retVal != 0) {
        return 1;
    }

    // Read mRNA fasta file
    retVal = mMRNASeqs.read_fasta(mMRNAFasta);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

} // namespace mikan
