#ifndef MK_INPUT_HPP_
#define MK_INPUT_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE
#include "mk_option.hpp"         // MKOptions
#include "mk_sequence.hpp"       // MKSequences

namespace mikan {

//
// Input data for mikan score
//
class MKInput {
public:
    // Define variables
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;

    // Define methods
    MKInput() {}

    mikan::TCharSet const &get_mirna_ids() { return mMiRNASeqs.get_ids(); }

    mikan::TRNASet const &get_mirna_seqs() { return mMiRNASeqs.get_seqs(); }

    mikan::TCharSet const &get_mrna_ids() { return mMRNASeqs.get_ids(); }

    mikan::TRNASet const &get_mrna_seqs() { return mMRNASeqs.get_seqs(); }

    // Method prototypes
    int load_seq_from_file();

    void set_file_names(seqan::CharString &pMiRNAFasta, seqan::CharString &pMRNAFasta);

    void set_options(MKOptions &opt);

private:
    // Define variables
    mikan::MKSequences mMiRNASeqs;
    mikan::MKSequences mMRNASeqs;
};

} // namespace mikan

#endif /* MK_INPUT_HPP_ */
