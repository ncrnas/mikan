#ifndef MK_TYPEDEF_HPP_
#define MK_TYPEDEF_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_const.hpp"

namespace mikan {

typedef seqan::RnaString TRNATYPE;
//typedef seqan::Rna5String TRNATYPE;

typedef seqan::StringSet<seqan::CharString> TCharSet;
typedef seqan::StringSet<TRNATYPE> TRNASet;

typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
typedef seqan::Finder<TIndexQGram> TFinder;

} // namespace mikan

#endif /* MK_TYPEDEF_HPP_ */
