#ifndef MK_SEED_SEQ_HPP_
#define MK_SEED_SEQ_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_const.hpp"
#include "mk_option.hpp"          // MKOptions

namespace mikan {

//
// Generate miRNA seeds
//
class MKSeedSeqs {
public:
    // Define variable
    seqan::String<bool> mEffectiveSeeds;
    mikan::TCharSet mSeedTypeDef;

    // Define methods
    MKSeedSeqs(mikan::MKOptions const &opts) : mOpts(opts) {
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
        mMMGULab = "MMGU";
        mMMLab = "MM";

        mFilterRedundant = true;
        mTSSVMMismatch = false;

        init_from_args();
        set_flags();
    }

    mikan::TRNAStr const &get_seed_seq(int pIdx) const { return mSeedSeqs[pIdx]; }

    mikan::TCharStr const &get_seed_type(int pIdx) const { return mSeedTypes[pIdx]; }

    unsigned get_mismatched_pos(int pIdx) { return mMisMatchPos[pIdx]; }

    // Method prototypes
    void init_from_args();

    void set_seed_type_def(mikan::TCharSet &pSeedTypeDef) { mSeedTypeDef = pSeedTypeDef; }

    void set_flags();

    void clear_seeds();

    int create_seed_seqs(mikan::TRNAStr const &pSeq);

    mikan::TRNAStr const get_mirna_seq() const { return mMiRNASeq; };

    void print_all();

protected:
    // Define variables
    mikan::MKOptions const &mOpts;
    mikan::TRNASet mSeedSeqs;
    mikan::TCharSet mSeedTypes;
    mikan::TSitePosSet mMisMatchPos;
    mikan::TRNAStr mMiRNASeq;
    mikan::TRNAStr mRNAChar;
    bool mSingleGU;
    bool mMultiGU;
    bool mMisMatch;
    bool mGUMisMatch;
    bool mBT;
    bool mBTM8;
    bool mBM;
    bool mLP;
    bool mOther;
    bool mAddInReverse;
    bool mFilterRedundant;
    bool mTSSVMMismatch;

    mikan::TCharStr mNMerLab;
    mikan::TCharStr mGUTLab;
    mikan::TCharStr mGUMLab;
    mikan::TCharStr mMultiGUTLab;
    mikan::TCharStr mMultiGUMLab;
    mikan::TCharStr mMMGULab;
    mikan::TCharStr mMMLab;

    unsigned nNumNewSeq;
    mikan::TRNASet mTmpSeedSeqs;
    mikan::TCharSet mTmpSeedTypes;
    mikan::TSitePosSet mTmpMisMatchPos;

    // Method prototypes
    int create_single_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_multi_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq, bool pIsGUMM = false, int pGUPos = 0);

    int create_gu_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_bt_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_bt_m8_seed_seqs();

    int create_bm_seed_seqs();

    virtual int create_other_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int check_redundant_seeds();

    int add_seed_seqs();

    void init_temp(unsigned pVecSize);

};

} // namespace mikan

#endif /* MK_SEED_SEQ_HPP_ */
