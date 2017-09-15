#ifndef MR3_ALIGN_HPP_
#define MR3_ALIGN_HPP_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "dp_core.hpp"            // MR3DPCore
#include "dp_score.hpp"           // MR3DPScore

namespace mr3as {

//
// Store alignments
//
class MR3Align {
public:
    // Define variables
    seqan::String<bool> mEffectiveSites;
    mikan::TCharSet mAlignMRNA;
    mikan::TCharSet mAlignBars;
    mikan::TCharSet mAlignMiRNA;

public:
    // Define methods
    MR3Align() {}

    int get_align_score(int pIdx) { return mAlignScores[pIdx]; }

    // Method prototype
    void clear_align();

    void resize_align(unsigned pSize);

    void align_seed(int pIdx, mikan::TRNAStr &pIMiRNASeedSeq, mikan::TRNAStr &pIMRNASeedSeq, int pMMpos);

    void align_3p(int pIdx, mikan::TRNAStr &pIMiRNA3pSeq, mikan::TRNAStr &pIMRNA3pSeq);

    void combine_alignments(int pIdx, mikan::TRNAStr const &pMiRNASeq, mikan::TRNAStr const &pMRNASeq, bool noA1);

    void get_mrna_seq(int pIdx, seqan::CharString &pStrMRNA);

    void init_3p_align(int pIdx);

private:
    seqan::String<int> mAlignSeedScores;
    seqan::String<int> mAlign3PScores;
    seqan::String<int> mAlignScores;

    mikan::TCharSet mAlignSeedMiRNA;
    mikan::TCharSet mAlignSeedMRNA;
    mikan::TCharSet mAlign3pMiRNA;
    mikan::TCharSet mAlign3pMRNA;
    seqan::String<int> mGapCountMiRNA;
    seqan::String<int> mGapCountMRNA;

private:
    void set_align_bars(int pIdx);

    mr3dp::MR3DPCore mDPCore;
    mr3dp::MR3DPScore mDPScore;

};

} // namespace mr3as

#endif /* MR3_ALIGN_HPP_ */
