#ifndef TSSVM_ALIGN_HPP_
#define TSSVM_ALIGN_HPP_

#include <tssvm_seed_site.hpp>      // TSSVMSeedSites
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>

// Extend SeqAn by a user-define scoring matrix.
namespace seqan {

// We have to create a new specialization of the ScoringMatrix_ class
// for amino acids.  For this, we first create a new tag.
struct ERNAGUWOB {
};

template<>
struct ScoringMatrixData_<int, Rna, ERNAGUWOB> {
    enum {
        VALUE_SIZE = ValueSize<Rna>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const *getData() {
        // The user defined data table.
        //      A   C   G   U
        // 	A  -3  -3  -3   5
        // 	C  -3  -3   5  -3
        // 	G  -3   5  -3   2
        // 	U   5  -3   2  -3
        static int const _data[TAB_SIZE] = {
                -3, -3, -3, 5,
                -3, -3, 5, -3,
                -3, 5, -3, 2,
                5, -3, 2, -3,
        };
        return _data;
    }
};

template<>
struct ScoringMatrixData_<int, Rna5, ERNAGUWOB> {
    enum {
        VALUE_SIZE = ValueSize<Rna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const *getData() {
        // The user defined data table.
        static int const _data[TAB_SIZE] = {
                -3, -3, -3, 5, -3,
                -3, -3, 5, -3, -3,
                -3, 5, -3, 2, -3,
                5, -3, 2, -3, -3,
                -3, -3, -3, -3, -3
        };
        return _data;
    }
};

}  // namespace seqan

namespace tssvm {

//
// Store miRNA and mRNA sequences and ids
//
template<class TRNAString>
class TSAlign {
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Constant values
    static const unsigned ALIGN_SEQ_LEN = 20;
    static const int GAP_OPEN_SCORE = -10;
    static const int GAP_EXTENT_SCORE = -7;

public:
    // Define methods
    TSAlign() : mScoreMatrix(GAP_EXTENT_SCORE, GAP_OPEN_SCORE) {}

    TCharSet const &get_align_bars() const { return mAlignBars; }

    TCharSet const &get_align_mrna() const { return mAlignMRNA; }

    TCharSet const &get_align_mirna() const { return mAlignMiRNA; }

    // Method prototypes
    int align_seq(TSSVMSeedSites<TRNAString> &pSeedSites, TRNAString const &pMiRNASeq,
                  TRNASet const &pMRNASeqs);

    void clear_alignments();

    void write_alignment(int pIdx);

private:
    TCharSet mAlignBars;
    TCharSet mAlignMRNA;
    TCharSet mAlignMiRNA;
    seqan::String<int> mAlignScores;

    typedef seqan::ScoreMatrix<typename seqan::Value<TRNAString>::Type, seqan::ERNAGUWOB> TScoreMat;
    seqan::Score<int, TScoreMat> mScoreMatrix;

    typedef seqan::Align<TRNAString, seqan::ArrayGaps> TAlign;
    typedef seqan::Gaps<TRNAString, seqan::ArrayGaps> TGap;

private:
    int set_mirna_seq_for_align(const seqan::CharString &pSeedType, TRNAString const &pMiRNASeq,
                                TRNAString &pMiRNAAlignSeq);

    int set_mrna_seq_for_align(const seqan::CharString &pSeedType, unsigned pSitePos, unsigned pMiRNALen,
                               const TRNAString &pMRNASeq, TRNAString &pMRNAAlignSeq);

    int set_addtional_sequences(TRNAString &pMiRNAAlignSeq, TRNAString &pMRNAAlignSeq);

    unsigned get_align_len(TAlign &pAlign);

    int set_align_mrna(TAlign &pAlign, unsigned pAlignLen, const seqan::CharString &pSeedType,
                       unsigned pMisMatchPos, unsigned pSitePos, const TRNAString &pMRNASeq, int pIdx);

    int set_align_mirna(TAlign &pAlign, unsigned pAlignLen, const seqan::CharString &pSeedType,
                        unsigned pMisMatchPos, const TRNAString &pMiRNASeq, int pIdx);

    int set_align_bars(int pIdx, unsigned pAlignLen);
};

} // namespace tssvm

#endif /* TSSVM_ALIGN_HPP_ */
