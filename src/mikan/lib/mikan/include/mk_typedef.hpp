#ifndef MK_TYPEDEF_HPP_
#define MK_TYPEDEF_HPP_

#include <seqan/sequence.h>

namespace mikan {

typedef seqan::RnaString TRNATYPE;
//typedef seqan::Rna5String TRNATYPE;

typedef seqan::StringSet<seqan::CharString> TCharSet;
typedef seqan::StringSet<TRNATYPE> TRNASet;

} // namespace mikan

#endif /* MK_TYPEDEF_HPP_ */
