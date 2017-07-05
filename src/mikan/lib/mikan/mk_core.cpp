#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif
#include <seqan/arg_parse.h>
#include <mk_inst_template.hpp>  // TRNATYPE
#include <mk_core.hpp>           // MR3CoreInput, MR3Core

#include <mr3_core.hpp>          // MR3CoreMain
#include <pita_core.hpp>          // PITACoreMain
#include <rh2_core.hpp>          // RH2CoreMain
#include <tm1_core.hpp>          // TM1CoreMain
#include <ts5_core.hpp>          // TS5CoreMain
#include <tssvm_core.hpp>           // TSSVMCoreMain

namespace mikan {

//
// MKCoreInput methods
//
template <class TRNAString>
void MKCoreInput<TRNAString>::init_from_args(MKOptions& opts)
{
    mMiRNAFasta = opts.mMiRNAFasta;
    mMRNAFasta = opts.mMRNAFasta;
}

template <class TRNAString>
int MKCoreInput<TRNAString>::load_seq_from_file()
{
    int retVal;

    // Read miRNA fasta file
    retVal = mMiRNASeqs.read_fasta(mMiRNAFasta);
    if (retVal != 0)
    {
        return 1;
    }

    // Read mRNA fasta file
    retVal = mMRNASeqs.read_fasta(mMRNAFasta);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}

// Explicit template instantiation
template class MKCoreInput<TRNATYPE>;

} // namespace mikan
