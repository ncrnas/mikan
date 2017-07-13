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
    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    MKSeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            mMRNASeqs(pMRNASeqs), mRNAIdx(pRNAIdx), mFinder(pFinder) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    mikan::TMRNAPosSet const &get_mrna_pos() const { return mMRNAPos; }

    mikan::TSitePosSet const &get_site_pos() const { return mSitePos; }

    mikan::TCharSet const &get_seed_types() const { return mSeedTypes; }

    mikan::TSitePosSet const &get_mismatched_pos() const { return mMisMatchPos; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::MKSeedSeqs &seedSeqs);

    void clear_pos();

protected:
    mikan::TMRNAPosSet mMRNAPos;
    mikan::TSitePosSet mSitePos;

    mikan::TRNASet const &mMRNASeqs;
    mikan::TCharSet mSeedTypes;
    mikan::TSitePosSet mMisMatchPos;

    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;

protected:
    virtual bool check_position(unsigned pMRNAPos, unsigned pSitePos);

    virtual void set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos, mikan::TRNAStr &pSeedSeq,
                                   seqan::CharString &pSeedType, unsigned pMisMatchPos, bool pEffectiveSite);

};

} // namespace mikan

#endif /* MK_SEED_SITE_HPP_ */
