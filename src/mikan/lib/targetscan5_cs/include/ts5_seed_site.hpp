#ifndef TS5_SEED_SITE_HPP_
#define TS5_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder

namespace ts5cs {

//
// Generate miRNA seeds
//
class TS5SeedSeqs {
public:
    // Define methods
    TS5SeedSeqs() {}

    mikan::TRNATYPE const &get_seed_seq() const { return mSeedSeqs[0]; }

    // Method prototypes
    int create_seed_seqs();

    void set_mirna_seq(mikan::TRNATYPE pSeq);

private:
    mikan::TRNASet mSeedSeqs;
    mikan::TRNATYPE mMiRNASeq;
};

//
// miRNA seed sites
//
class TS5SeedSites {
public:
    // Define methods
    TS5SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) : 
            mRNAIdx(pRNAIdx), mFinder(pFinder) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::TRNATYPE const &pMiRNA);

    void clear_pos();

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;
};

} // namespace ts5cs

#endif /* TS5_SEED_SITE_HPP_ */
