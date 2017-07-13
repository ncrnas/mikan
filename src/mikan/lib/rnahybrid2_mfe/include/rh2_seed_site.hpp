#ifndef RH2_SEED_SITE_HPP_
#define RH2_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"        // MKSeedSeqs

namespace rh2mfe {

//
// Generate miRNA seeds
//
class RH2SeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    RH2SeedSeqs() : MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet &pSeedTypeDef);
};

//
// miRNA seed sites
//
class RH2SeedSites {
public:
    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 1;
    static const unsigned MIN_DIST_UTR_END = 0;
    static const unsigned SEED_LEN = 6;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    RH2SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    seqan::StringSet<seqan::CharString> const &get_seed_types() const { return mSeedTypes; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::TRNAStr const &pMiRNA, mikan::TCharSet &pSeedDef);

    void clear_pos();

    void set_new_seed_type(seqan::CharString &pCurSeedType, seqan::CharString &pSeedDef,
                           unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr const &pMiRNA, bool &pEffectiveSite);

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;
    mikan::TRNASet const &mMRNASeqs;

};

} // namespace rh2mfe

#endif /* RH2_SEED_SITE_HPP_ */
