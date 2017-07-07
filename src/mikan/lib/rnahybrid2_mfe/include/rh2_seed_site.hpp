#ifndef RH2_SEED_SITE_HPP_
#define RH2_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>

namespace rh2mfe {

//
// Generate miRNA seeds
//
template<class TRNAString>
class RH2SeedSeqs {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Constant values
    static const unsigned SEED_LEN = 6;

    // Define variables
    seqan::String<bool> mEffectiveSeeds;

public:
    // Define methods
    RH2SeedSeqs() {}

    TRNAString const &get_seed_seq(int i) const { return mSeedSeqs[i]; }

    seqan::CharString const &get_seed_type(int i) const { return mSeedTypes[i]; }

    unsigned get_mismatched_pos(int) { return 0; }

    // Method prototypes
    int create_seed_seqs(seqan::CharString &pSeedType, seqan::CharString &pOverlapDef);

    void set_mirna_seq(TRNAString pSeq);

private:
    TRNASet mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    TRNAString mMiRNASeq;

private:
    int create_nmer_seed_seqs(TRNAString &pSeedSeq, seqan::CharString &pSeedDef);

    int create_single_guwobble_seed_seqs(TRNAString &pSeedSeq, seqan::CharString &pGUT, seqan::CharString &pGUM);

    int create_multi_guwobble_seed_seqs(TRNAString &pSeedSeq, seqan::CharString &pGUT, seqan::CharString &pGUM);

    int check_redundant_seeds(seqan::CharString &pOverlapDef);

};

//
// miRNA seed sites
//
template<class TRNAString>
class RH2SeedSites {
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<6> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 1;
    static const unsigned MIN_DIST_UTR_END = 0;
    static const unsigned SEED_LEN = 6;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    RH2SeedSites(TIndexQGram &pRNAIdx, TFinder &pFinder, TRNASet const &pMRNASeqs) :
            mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    seqan::StringSet<seqan::CharString> const &get_seed_types() const { return mSeedTypes; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(TRNAString const &pMiRNA, seqan::CharString &pSeedDef,
                        seqan::CharString &pOverlapDef);

    void clear_pos();

    void set_new_seed_type(seqan::CharString &pCurSeedType, seqan::CharString &pSeedDef,
                           unsigned pMRNAPos, unsigned pSitePos, TRNAString const &pMiRNA, bool &pEffectiveSite);

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    TIndexQGram &mRNAIdx;
    TFinder &mFinder;
    TRNASet const &mMRNASeqs;

};

} // namespace rh2mfe

#endif /* RH2_SEED_SITE_HPP_ */
