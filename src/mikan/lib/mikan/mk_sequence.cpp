#include "mk_typedef.hpp"  // TRNATYPE
#include "mk_sequence.hpp"       // MKSequences
#include <seqan/seq_io.h>

using namespace seqan;

namespace mikan {

//
// MKSequences methods
//
template<class TRNAString>
int MKSequences<TRNAString>::read_fasta(CharString const &pFasta) {
    CharString id;
    CharString seq;

    SequenceStream seqStream(toCString(pFasta));
    if (!isGood(seqStream)) {
        std::cerr << "ERROR: Could not open the file!" << std::endl;
        return 1;
    }

    while (!atEnd(seqStream)) {
        if (readRecord(id, seq, seqStream) != 0) {
            std::cerr << "ERROR: Could not read from " << toCString(pFasta) << "!" << std::endl;
            return 1;
        }

        toUpper(seq);
        for (unsigned i = 0; i < length(seq); ++i) {
            if (seq[i] == 'T') {
                seq[i] = 'U';
            }
        }
        if (length(seq) > 0) {
            appendValue(mSeqIds, id);
            appendValue(mSeqs, seq);
        }

        if (mMaxLen < (int) length(seq)) {
            mMaxLen = (int) length(seq);
        }
    }

    return 0;
}

// Explicit template instantiation
template
class MKSequences<TRNATYPE>;

} // namespace mikan
