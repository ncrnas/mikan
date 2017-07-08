#ifndef TM1_ALIGN_HPP_
#define TM1_ALIGN_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_align.hpp"         // TM1Alignment
#include "tm1_seed_site.hpp"     // TM1SeedSites

namespace tm1p {

//
//  miRNA:mRNA alignment
//
class TM1Alignment {
public:
    // Define methods
    TM1Alignment() {}

    mikan::TCharSet const &get_align_bars() const { return mAlignBars; }

    mikan::TCharSet const &get_align_mrna() const { return mAlignMRNA; }

    mikan::TCharSet const &get_align_mirna() const { return mAlignMiRNA; }

    // Method prototypes
    void
    align_all(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs, TM1SeedSites const &pSeedSites,
              seqan::String<int> const &pA1Pos);

    void resize_alignments(unsigned pSize);

    void clear_alignments();

    void write_alignment(int pIdx) const;

private:
    mikan::TCharSet mAlignBars;
    mikan::TCharSet mAlignMRNA;
    mikan::TCharSet mAlignMiRNA;

};

} // namespace tm1p

#endif /* TM1_ALIGN_HPP_ */
