#ifndef TS5_ALIGN_HPP_
#define TS5_ALIGN_HPP_

#include <seqan/sequence.h>

namespace ts5cs {

//
//  miRNA:mRNA alignment
//
template<class TRNAString>
class TS5Alignment {
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

public:
    // Define methods
    TS5Alignment() {}

    TCharSet const &get_align_bars() const { return mAlignBars; }

    TCharSet const &get_align_mrna() const { return mAlignMRNA; }

    TCharSet const &get_align_mirna() const { return mAlignMiRNA; }

    // Method prototypes
    int align_seed(unsigned pMRNAIdx, const seqan::CharString &pSeedType, const TRNAString &pMiRNASeq,
                   const TRNAString &pMRNASeq, unsigned pSitePos);

    int align_3p_part(unsigned pMRNAIdx, const seqan::CharString &pSeedType, const TRNAString &pMiRNAThreePrime,
                      TRNAString &pMRNAThreePrime, seqan::String<int> &pMatchLen, seqan::String<int> &pMiRNAPos,
                      seqan::String<int> &pMRNAPos, float pScore, unsigned pMatchedIdx);

    void resize_alignments(unsigned pSize);

    void clear_alignments();

    void write_alignment(int pIdx) const;

private:
    TCharSet mAlignBars;
    TCharSet mAlignMRNA;
    TCharSet mAlignMiRNA;

private:
    int align_no_3p_part(unsigned pMRNAIdx, TRNAString &pMRNAThreePrime, TRNAString &pMiRNAThreePrime,
                         unsigned pAlignLen);

    void set_alignment(unsigned pMRNAIdx, TRNAString &pSeqThreePrime, seqan::String<int> &pSeqPos1,
                       seqan::String<int> &pSeqPos2, unsigned pMatchedIdx, TCharSet &pAlignSeq);

    void set_align_bars(unsigned pMRNAIdx, seqan::String<int> &pMatchLen, seqan::String<int> &pMRNAPos,
                        seqan::String<int> &pMiRNAPos, unsigned pMatchedIdx);

    void set_mismatch(unsigned pMRNAIdx, TCharSet &pAlignSeq);

    void trim_mismatch(unsigned pMRNAIdx);

};

} // namespace ts5cs

#endif /* TS5_ALIGN_HPP_ */
