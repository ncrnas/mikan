#ifndef MK_SEQUENCE_HPP_
#define MK_SEQUENCE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>

namespace mikan {

//
// Store miRNA and mRNA sequences and ids
//
template<class TRNAString>
class MKSequences {
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Define methods
    MKSequences() : mMaxLen(0) {}

    unsigned get_length() const { return length(mSeqIds); }

    seqan::CharString const &get_seq_id(unsigned const pIdx) const { return mSeqIds[pIdx]; }

    TRNAString const &get_seq(unsigned const pIdx) const { return mSeqs[pIdx]; }

    TCharSet const &get_ids() const { return mSeqIds; }

    TRNASet const &get_seqs() const { return mSeqs; }

    int get_max_seq_len() { return mMaxLen; }

    // Method prototypes
    int read_fasta(seqan::CharString const &pFasta);

private:
    TCharSet mSeqIds;
    TRNASet mSeqs;
    int mMaxLen;
};

} // namespace mikan

#endif /* MK_SEQUENCE_HPP_ */
