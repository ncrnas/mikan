#ifndef TM1_SEED_SITE_HPP_
#define TM1_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder

namespace tm1p {

//
// Generate miRNA seeds
//
class TM1SeedSeqs {
public:
    // Define variables
    seqan::String<bool> mEffectiveSeeds;

public:
    // Define methods
    TM1SeedSeqs() {}

    mikan::TRNATYPE const &get_seed_seq(int i) const { return mSeedSeqs[i]; }

    seqan::CharString const &get_seed_type(int i) const { return mSeedTypes[i]; }

    // Method prototypes
    int create_seed_seqs();

    void set_mirna_seq(mikan::TRNATYPE pSeq);

private:
    mikan::TRNASet mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    mikan::TRNATYPE mMiRNASeq;

private:
    int create_nmer_seed_seqs(mikan::TRNATYPE &pSeedSeq);

    int create_single_guwobble_seed_seqs(mikan::TRNATYPE &pSeedSeq);

    int create_lp_seed_seqs(mikan::TRNATYPE &pSeedSeq, seqan::CharString &pSeedType);

    int check_redundant_seeds();
};

//
// miRNA seed sites
//
class TM1SeedSites {
public:
    // Define types
    typedef seqan::StringSet<mikan::TRNATYPE> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<6> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 1;
    static const unsigned MIN_DIST_UTR_END = 0;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    TM1SeedSites(TIndexQGram &pRNAIdx, TFinder &pFinder, TRNASet const &pMRNASeqs) :
            mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    seqan::StringSet<seqan::CharString> const &get_seed_types() const { return mSeedTypes; }

    bool is_m8_match_gu(int i) { return (mM8Match[i] || mM8GU[i]); }

    bool is_m8_match(int i) { return (mM8Match[i]); }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::TRNATYPE const &pMiRNA);

    void clear_pos();

    int get_seed_len(int pIdx);

    int get_seed_start_pos(int pIdx);

    int get_seed_start_pos2(int pIdx);

    int get_seed_end_pos(int pIdx);

    int get_seed_end_pos2(int pIdx);

    int get_length_to_cds(int pIdx);

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;
    mikan::TRNASet const &mMRNASeqs;

    seqan::String<bool> mM8Match;
    seqan::String<bool> mM8GU;
    seqan::String<bool> mM1A;
    seqan::String<bool> mM1Match;
    seqan::String<bool> mM1GU;
    seqan::String<unsigned> mMRNASeqLen;
    seqan::String<unsigned> mM8Pos;

private:
    void set_new_seed_type(seqan::CharString &pCurSeedType, unsigned pMRNAPos, unsigned pSitePos,
                           mikan::TRNATYPE const &pMiRNA, bool &pEffectiveSite);

    void get_mx_match(mikan::TRNATYPE const &pMiRNASeq, mikan::TRNATYPE const &pMiRNACompSeq, mikan::TRNATYPE const &pMRNASeq,
                      unsigned pSitePos, int pMx, bool &pMatch, bool &pGU, bool &pNoMx, bool &pIsA);

    void get_match_count(mikan::TRNATYPE const &pMiRNASeq, mikan::TRNATYPE const &pMiRNACompSeq, mikan::TRNATYPE const &pMRNASeq,
                         unsigned pSitePos, int pMx1, int pMx2, int &pMatchCount, int &pGUCount);
};


} // namespace tm1p

#endif /* TM1_SEED_SITE_HPP_ */
