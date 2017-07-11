#ifndef MR3_SEED_SITE_HPP_
#define MR3_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seqs.hpp"       // MKSeedSeqs

namespace mr3as {

//
// Generate miRNA seeds
//
class MR3SeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    MR3SeedSeqs() : MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet pSeedTypeDef);
};

//
// miRNA seed sites
//
class MR3SeedSites {
public:
    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 7;
    static const unsigned MIN_DIST_UTR_END = 0;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const bool FORCE_LAST_MATCH = false;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    MR3SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    seqan::StringSet<seqan::CharString> const &get_seed_types() const { return mSeedTypes; }

    seqan::String<int> const &get_mismatched_pos() const { return mMisMatchPos; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::TRNAStr const &pMiRNA, mikan::TCharSet &pSeedDef);

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
                                 bool pMatchMx8, bool pMatchMx9, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_single_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                 int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                 bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos,
                                 seqan::CharString &pNewSeedType);

    void set_multiple_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                   bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void set_gu_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                                   int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
                                   bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos,
                                   seqan::CharString &pNewSeedType);

    void set_bt_seed_type(unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, unsigned pMisMatchPos,
                          seqan::CharString &pNewSeedType);

    void set_6mer_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
                            bool pMatchMx8, bool pMatchMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType);

};

} // namespace mr3as

#endif /* MR3_SEED_SITE_HPP_ */
