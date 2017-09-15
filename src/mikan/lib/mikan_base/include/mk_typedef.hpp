#ifndef MK_TYPEDEF_HPP_
#define MK_TYPEDEF_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_const.hpp"

namespace mikan {

typedef seqan::CharString TCharStr;
typedef seqan::RnaString TRNAStr;
//typedef seqan::Rna5String TRNAStr;

typedef seqan::StringSet<mikan::TCharStr> TCharSet;
typedef seqan::StringSet<TRNAStr> TRNASet;

typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
typedef seqan::Finder<TIndexQGram> TFinder;

typedef seqan::String<unsigned> TMRNAPosSet;
typedef seqan::String<unsigned> TSitePosSet;
typedef seqan::String<int> TMismatchSet;
typedef seqan::String<float> TScoreSet;

} // namespace mikan

#endif /* MK_TYPEDEF_HPP_ */
