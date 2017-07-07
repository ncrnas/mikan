#ifndef MK_INPUT_HPP_
#define MK_INPUT_HPP_

#include "mk_typedef.hpp"  // TRNATYPE
#include "mk_option.hpp"         // MKOptions
#include "mk_sequence.hpp"       // MKSequences
#include <seqan/sequence.h>

namespace mikan {

//
// Input data for mikan score
//
class MKInput {
public:

    // Declare variables
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;

public:
    // Define methods
    MKInput() {}

    TCharSet const &get_mirna_ids() { return mMiRNASeqs.get_ids(); }

    TRNASet const &get_mirna_seqs() { return mMiRNASeqs.get_seqs(); }

    TCharSet const &get_mrna_ids() { return mMRNASeqs.get_ids(); }

    TRNASet const &get_mrna_seqs() { return mMRNASeqs.get_seqs(); }

    // Method prototypes
    int load_seq_from_file();

    void set_file_names(seqan::CharString &pMiRNAFasta, seqan::CharString &pMRNAFasta);

    void set_options(MKOptions &opt);

private:
    MKSequences mMiRNASeqs;
    MKSequences mMRNASeqs;
};

} // namespace mikan

#endif /* MK_INPUT_HPP_ */
