#ifndef TS5_SEED_SITE_HPP_
#define TS5_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>

namespace ts5cs {

//
// Generate miRNA seeds
//
template<class TRNAString>
class TS5SeedSeqs {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

public:
    // Define methods
    TS5SeedSeqs() {}

    TRNAString const &get_seed_seq() const { return mSeedSeqs[0]; }

    // Method prototypes
    int create_seed_seqs();

    void set_mirna_seq(TRNAString pSeq);

private:
    TRNASet mSeedSeqs;
    TRNAString mMiRNASeq;
};

//
// miRNA seed sites
//
template<class TRNAString, int SEEDLEN = 6>
class TS5SeedSites {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<SEEDLEN> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

public:
    // Define methods
    TS5SeedSites(TIndexQGram &pRNAIdx, TFinder &pFinder) : mRNAIdx(pRNAIdx), mFinder(pFinder) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(TRNAString const &pMiRNA);

    void clear_pos();

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    TIndexQGram &mRNAIdx;
    TFinder &mFinder;
};

} // namespace ts5cs

#endif /* TS5_SEED_SITE_HPP_ */
