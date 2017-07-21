#ifndef MK_SEQUENCE_HPP_
#define MK_SEQUENCE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"

namespace mikan {

//
// Store miRNA and mRNA sequences and ids
//
class MKSequences {
public:
    // Define methods
    explicit MKSequences() : mMaxLen(0) {}

    unsigned get_length() const { return length(mSeqIds); }

    seqan::CharString const &get_seq_id(unsigned const pIdx) const { return mSeqIds[pIdx]; }

    mikan::TRNAStr const &get_seq(unsigned const pIdx) const { return mSeqs[pIdx]; }

    mikan::TCharSet const &get_ids() const { return mSeqIds; }

    mikan::TRNASet const &get_seqs() const { return mSeqs; }

    int get_max_seq_len() { return mMaxLen; }

    // Method prototype
    int read_fasta(seqan::CharString const &pFasta);

private:
    // Define variables
    mikan::TCharSet mSeqIds;
    mikan::TRNASet mSeqs;
    int mMaxLen;

};

} // namespace mikan

#endif /* MK_SEQUENCE_HPP_ */
