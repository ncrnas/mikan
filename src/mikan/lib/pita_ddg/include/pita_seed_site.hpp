#ifndef PITA_SEED_SITE_HPP_
#define PITA_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder

namespace ptddg {

//
// Generate miRNA seeds
//
class PITASeedSeqs {
public:

    // Define variables
    seqan::String<bool> mEffectiveSeeds;

public:
    // Define methods
    PITASeedSeqs() {}

    mikan::TRNAStr const &get_seed_seq(int i) const { return mSeedSeqs[i]; }

    seqan::CharString const &get_seed_type(int i) const { return mSeedTypes[i]; }

    unsigned get_mismatched_pos(int i) { return mMisMatchPos[i]; }

    // Method prototypes
    int create_seed_seqs(seqan::StringSet<seqan::CharString> &pSeedType);

    void set_mirna_seq(mikan::TRNAStr pSeq);

private:
    mikan::TRNASet mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    seqan::String<unsigned> mMisMatchPos;
    mikan::TRNAStr mMiRNASeq;

private:
    int create_nmer_seed_seqs(mikan::TRNAStr &pSeedSeq, seqan::StringSet<seqan::CharString> &pSeedDef);

    int create_single_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_multi_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq, bool pIsGUMM = false, int pGUPos = 0);

    int create_gu_mismatch_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int check_redundant_seeds();

};

//
// miRNA seed sites
//
class PITASeedSites {
public:

    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 22;
    static const unsigned MIN_DIST_UTR_END = 1;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const bool FORCE_LAST_MATCH = true;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    PITASeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    seqan::StringSet<seqan::CharString> const &get_seed_types() const { return mSeedTypes; }

    seqan::String<int> const &get_mismatched_pos() const { return mMisMatchPos; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::TRNAStr const &pMiRNA, seqan::StringSet<seqan::CharString> &pSeedDef);

    void clear_pos();

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    seqan::String<int> mMisMatchPos;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;
    mikan::TRNASet const &mMRNASeqs;

private:
    void set_new_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                           unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, unsigned pMisMatchPos,
                           bool &pEffectiveSite);

    void set_mx_matches(unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, int pMx,
                        bool &pMatchMx, bool &pGutMx, bool &pGumMx);

    void set_stringent_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 bool pMatchMx1, bool pMatchMx2, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_single_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                 bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_multiple_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                   bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void set_gu_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx1, bool pMatchMx2, bool pGutMx1, bool pGutMx2,
                                   bool pGumMx1, bool pGumMx2, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_6mer_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                            bool pMatchMx1, bool pMatchMx2, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType);

};

} // namespace ptddg

#endif /* PITA_SEED_SITE_HPP_ */
