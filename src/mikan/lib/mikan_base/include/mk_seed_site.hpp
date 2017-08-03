#ifndef MK_SEED_SITE_HPP_
#define MK_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"       // MKSeedSeqs

namespace mikan {

//
// miRNA seed sites
//
class MKSeedSites {
public:
    // Constant value
    static const unsigned INDEXED_SEQ_LEN = SEEDLEN;

    // Define variable
    seqan::String<bool> mEffectiveSites;

    // Define methods
    explicit MKSeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            mMRNASeqs(pMRNASeqs), mRNAIdx(pRNAIdx), mFinder(pFinder) {

        mMinToCDS = 0;
        mMinToEnd = 0;

        mUpdatePos = true;
    }

    unsigned get_length() const { return seqan::length(mSitePos); }

    mikan::TMRNAPosSet const &get_mrna_pos() const { return mMRNAPos; }

    mikan::TSitePosSet const &get_site_pos() const { return mSitePos; }

    mikan::TCharSet const &get_seed_types() const { return mSeedTypes; }

    mikan::TMismatchSet const &get_mismatched_pos() const { return mMisMatchPos; }

    mikan::TSitePosSet const &get_site_pos_s1() const { return mS1Pos; }

    mikan::TSitePosSet const &get_site_pos_s8() const { return mS8Pos; }

    unsigned get_mrna_seq_len(int pIdx) { return length(mMRNASeqs[mMRNAPos[pIdx]]); }

    virtual int get_seed_len(int) { return INDEXED_SEQ_LEN; }

    virtual int get_seed_start(int pIdx) { return mSitePos[pIdx]; }

    virtual int get_seed_end(int pIdx) { return mSitePos[pIdx] + INDEXED_SEQ_LEN + 1; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::MKSeedSeqs &seedSeqs, mikan::TCharSet &pSeedTypeDef);

    virtual void clear_pos();

    void print_all();

protected:
    // Define variables
    mikan::TMRNAPosSet mMRNAPos;
    mikan::TSitePosSet mSitePos;
    mikan::TCharSet mSeedTypes;
    mikan::TMismatchSet mMisMatchPos;

    mikan::TRNASet const &mMRNASeqs;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;

    unsigned mMinToCDS;
    unsigned mMinToEnd;

    bool mUpdatePos;
    mikan::TSitePosSet mS1Pos;
    mikan::TSitePosSet mS8Pos;

    // Define method
    virtual bool check_position_2(unsigned, unsigned, seqan::CharString &) { return true; }

    // Method prototypes
    virtual bool check_position_1(unsigned pMRNAPos, unsigned pSitePos, seqan::CharString &pSeedType);

    virtual bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                                   mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                                   seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

    void set_mx_matches(unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, int pMx,
                        bool &pNoMx, bool &pMatchMx, bool &pGutMx, bool &pGumMx, bool &pIsA);

};

} // namespace mikan

#endif /* MK_SEED_SITE_HPP_ */
