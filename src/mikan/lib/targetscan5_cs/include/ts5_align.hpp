#ifndef TS5_ALIGN_HPP_
#define TS5_ALIGN_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder

namespace ts5cs {

//
//  miRNA:mRNA alignment
//
class TS5Alignment {
public:
    // Define methods
    TS5Alignment() {}

    mikan::TCharSet const &get_align_bars() const { return mAlignBars; }

    mikan::TCharSet const &get_align_mrna() const { return mAlignMRNA; }

    mikan::TCharSet const &get_align_mirna() const { return mAlignMiRNA; }

    // Method prototypes
    int align_seed(unsigned pMRNAIdx, const seqan::CharString &pSeedType, const mikan::TRNAStr &pMiRNASeq,
                   const mikan::TRNAStr &pMRNASeq, unsigned pSitePos);

    int align_3p_part(unsigned pMRNAIdx, const seqan::CharString &pSeedType, const mikan::TRNAStr &pMiRNAThreePrime,
                      mikan::TRNAStr &pMRNAThreePrime, seqan::String<int> &pMatchLen, seqan::String<int> &pMiRNAPos,
                      seqan::String<int> &pMRNAPos, float pScore, unsigned pMatchedIdx);

    void resize_alignments(unsigned pSize);

    void clear_alignments();

    void write_alignment(int pIdx) const;

private:
    mikan::TCharSet mAlignBars;
    mikan::TCharSet mAlignMRNA;
    mikan::TCharSet mAlignMiRNA;

private:
    int align_no_3p_part(unsigned pMRNAIdx, mikan::TRNAStr &pMRNAThreePrime, mikan::TRNAStr &pMiRNAThreePrime,
                         unsigned pAlignLen);

    void set_alignment(unsigned pMRNAIdx, mikan::TRNAStr &pSeqThreePrime, seqan::String<int> &pSeqPos1,
                       seqan::String<int> &pSeqPos2, unsigned pMatchedIdx, mikan::TCharSet &pAlignSeq);

    void set_align_bars(unsigned pMRNAIdx, seqan::String<int> &pMatchLen, seqan::String<int> &pMRNAPos,
                        seqan::String<int> &pMiRNAPos, unsigned pMatchedIdx);

    void set_mismatch(unsigned pMRNAIdx, mikan::TCharSet &pAlignSeq);

    void trim_mismatch(unsigned pMRNAIdx);

};

} // namespace ts5cs

#endif /* TS5_ALIGN_HPP_ */
