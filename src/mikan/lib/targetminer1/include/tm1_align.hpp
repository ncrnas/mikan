#ifndef TM1_ALIGN_HPP_
#define TM1_ALIGN_HPP_

#include <tm1_align.hpp>         // TM1Alignment
#include <tm1_seed_site.hpp>     // TM1SeedSites
#include <seqan/sequence.h>

namespace tm1p {

//
//  miRNA:mRNA alignment
//
template<class TRNAString>
class TM1Alignment {
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::String<unsigned> TSitePos;

public:
    // Define methods
    TM1Alignment() {}

    TCharSet const &get_align_bars() const { return mAlignBars; }

    TCharSet const &get_align_mrna() const { return mAlignMRNA; }

    TCharSet const &get_align_mirna() const { return mAlignMiRNA; }

    // Method prototypes
    void
    align_all(TRNAString const &pMiRNASeq, TRNASet const &pMRNASeqs, TM1SeedSites <TRNAString> const &pSeedSites,
              seqan::String<int> const &pA1Pos);

    void resize_alignments(unsigned pSize);

    void clear_alignments();

    void write_alignment(int pIdx) const;

private:
    TCharSet mAlignBars;
    TCharSet mAlignMRNA;
    TCharSet mAlignMiRNA;

};

} // namespace tm1p

#endif /* TM1_ALIGN_HPP_ */
