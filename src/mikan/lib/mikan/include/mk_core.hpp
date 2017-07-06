#ifndef MK_CORE_HPP_
#define MK_CORE_HPP_

#include <mk_inst_template.hpp>  // TRNATYPE
#include <mk_option.hpp>         // MKOptions
#include <mk_sequence.hpp>       // MKSequences
#include <seqan/sequence.h>

namespace mikan {

//
// Input data for mikan score
//
template<class TRNAString>
class MKCoreInput {
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Declare variables
    seqan::CharString &mMiRNAFasta;
    seqan::CharString &mMRNAFasta;

public:
    // Define methods
    MKCoreInput(seqan::CharString &pMiRNAFasta, seqan::CharString &pMRNAFasta) :
            mMiRNAFasta(pMiRNAFasta), mMRNAFasta(pMRNAFasta){}

    TCharSet const &get_mirna_ids() { return mMiRNASeqs.get_ids(); }

    TRNASet const &get_mirna_seqs() { return mMiRNASeqs.get_seqs(); }

    TCharSet const &get_mrna_ids() { return mMRNASeqs.get_ids(); }

    TRNASet const &get_mrna_seqs() { return mMRNASeqs.get_seqs(); }

    // Method prototypes
    int load_seq_from_file();

private:
    MKSequences<TRNAString> mMiRNASeqs;
    MKSequences<TRNAString> mMRNASeqs;
};

} // namespace mikan

#endif /* MK_CORE_HPP_ */
