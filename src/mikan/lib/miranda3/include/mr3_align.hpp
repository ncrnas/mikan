#ifndef MR3_ALIGN_HPP_
#define MR3_ALIGN_HPP_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>

// Extend SeqAn by a user-define scoring matrix.
namespace seqan {

// We have to create a new specialization of the ScoringMatrix_ class
// for amino acids.  For this, we first create a new tag.
struct MR3RNA3P {};

template <>
struct ScoringMatrixData_<int, Rna5, MR3RNA3P> {
    enum {
        VALUE_SIZE = ValueSize<Rna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        // The user defined data table.
        //      A   C   G   U
        //  A  -3  -3  -3   5
        //  C  -3  -3   5  -3
        //  G  -3   5  -3   1
        //  U   5  -3   1  -3
        static int const _data[TAB_SIZE] = {
                -3,   -3,   -3,    5, -140,
                -3,   -3,    5,   -3, -140,
                -3,    5,   -3,    1, -140,
                 5,   -3,    1,   -3, -140,
              -140, -140, -140, -140,  140
        };
        return _data;
    }
};

}  // namespace seqan

namespace mr3as{

//
// Store alignments
//
template <class TRNAString>
class MR3Align
{
public:
    // Constant values
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const int GAP_OPEN_SCORE = -9;
    static const int GAP_EXTENT_SCORE = -4;
    static const int SEED_MISMATCH = -12;

    // Define variables
    seqan::String<bool> mEffectiveSites;
    seqan::StringSet<seqan::CharString> mAlignMRNA;
    seqan::StringSet<seqan::CharString> mAlignBars;
    seqan::StringSet<seqan::CharString> mAlignMiRNA;

public:
    // Define methods
    MR3Align(): mScoreMatrix3P(GAP_EXTENT_SCORE, GAP_OPEN_SCORE){}
    int get_align_score(int pIdx){return mAlignSeedScores[pIdx] + mAlign3PScores[pIdx];}

    // Method prototype
    void clear_align();
    void resize_align(unsigned pSize);
    void align_seed(int pIdx, TRNAString &pIMiRNASeedSeq, TRNAString &pIMRNASeedSeq, int pMMpos);
    void align_3p(int pIdx, seqan::Rna5String &pIMiRNA3pSeq, seqan::Rna5String &pIMRNA3pSeq);
    void combine_alignments(int pIdx, TRNAString const &pMiRNASeq, TRNAString const &pMRNASeq);
    void get_mrna_seq(int pIdx, TRNAString& pStrMRNA);
    void init_3p_align(int pIdx);

private:
    typedef seqan::Align<seqan::Rna5String, seqan::ArrayGaps> TAlign;
    typedef seqan::Gaps<seqan::Rna5String, seqan::ArrayGaps> TGap;
    typedef seqan::ScoreMatrix<typename seqan::Value<seqan::Rna5String>::Type, seqan::MR3RNA3P> TScore3PMat;

    seqan::Score<int, TScore3PMat> mScoreMatrix3P;

    TAlign mAign3P;
    seqan::String<int> mAlignSeedScores;
    seqan::String<int> mAlign3PScores;
    seqan::String<int> mAlignScores;

    seqan::StringSet<seqan::CharString> mAlignSeedMiRNA;
    seqan::StringSet<seqan::CharString> mAlignSeedMRNA;
    seqan::StringSet<seqan::CharString> mAlign3pMiRNA;
    seqan::StringSet<seqan::CharString> mAlign3pMRNA;
    seqan::String<int> mGapCount3pMiRNA;
    seqan::String<int> mGapCount3pMRNA;

private:
    void set_align_bars(int pIdx);
};

} // namespace mr3as

#endif /* MR3_ALIGN_HPP_ */
