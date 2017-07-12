#ifndef MK_SEED_SEQ_HPP_
#define MK_SEED_SEQ_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_const.hpp"

namespace mikan {

//
// Generate miRNA seeds
//
class MKSeedSeqs {
public:
    // Define variables
    seqan::String<bool> mEffectiveSeeds;

public:
    // Define methods
    MKSeedSeqs() {
        resize(mRNAChar, 4);
        mRNAChar[0] = 'A';
        mRNAChar[1] = 'C';
        mRNAChar[2] = 'G';
        mRNAChar[3] = 'U';

        init_temp(TEMP_SEED_SEQ_SIZE);
        mNMerLab = "6mer";
        mGUTLab = "GUT";
        mGUMLab = "GUM";
        mMultiGUTLab = "GU+";
        mMultiGUMLab = "GU+";
    }

    mikan::TRNAStr const &get_seed_seq(int i) const { return mSeedSeqs[i]; }

    seqan::CharString const &get_seed_type(int i) const { return mSeedTypes[i]; }

    unsigned get_mismatched_pos(int i) { return mMisMatchPos[i]; }

    // Method prototypes
    int create_seed_seqs();

    void set_flags(mikan::TCharSet &pSeedTypeDef);

    void set_mirna_seq(mikan::TRNAStr pSeq);

    void print_all();

protected:
    mikan::TRNASet mSeedSeqs;
    mikan::TCharSet mSeedTypes;
    seqan::String<unsigned> mMisMatchPos;
    mikan::TRNAStr mMiRNASeq;
    mikan::TRNAStr mRNAChar;
    bool mSingleGU;
    bool mMultiGU;
    bool mMisMatch;
    bool mGUMisMatch;
    bool mBT;
    bool mBM;
    bool mLP;
    bool mOther;
    bool mAddInReverse;

    seqan::CharString mNMerLab;
    seqan::CharString mGUTLab;
    seqan::CharString mGUMLab;
    seqan::CharString mMultiGUTLab;
    seqan::CharString mMultiGUMLab;

    unsigned nNumNewSeq;
    mikan::TRNASet mTmpSeedSeqs;
    mikan::TCharSet mTmpSeedTypes;
    seqan::String<unsigned> mTmpMisMatchPos;

protected:
    int create_single_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_multi_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq, bool pIsGUMM = false, int pGUPos = 0);

    int create_gu_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_bt_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_bm_seed_seqs();

    int create_lp_seed_seqs(mikan::TRNAStr &pSeedSeq);

    virtual int create_other_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int check_redundant_seeds();

    int add_seed_seqs();

    void init_temp(unsigned pVecSize);
};

} // namespace mikan

#endif /* MK_SEED_SEQ_HPP_ */
