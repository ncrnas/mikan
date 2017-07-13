#ifndef TM1_SEED_SITE_HPP_
#define TM1_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"       // MKSeedSeqs

namespace tm1p {

//
// Generate miRNA seeds
//
class TM1SeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    TM1SeedSeqs(): MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet &pSeedTypeDef);

protected:
    virtual int create_other_seed_seqs(mikan::TRNAStr &pSeedSeq);
};

//
// miRNA seed sites
//
class TM1SeedSites {
public:
    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 1;
    static const unsigned MIN_DIST_UTR_END = 0;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    TM1SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    seqan::StringSet<seqan::CharString> const &get_seed_types() const { return mSeedTypes; }

    bool is_m8_match_gu(int i) { return (mM8Match[i] || mM8GU[i]); }

    bool is_m8_match(int i) { return (mM8Match[i]); }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::TRNAStr const &pMiRNA);

    void clear_pos();

    int get_seed_len(int pIdx);

    int get_seed_start_pos(int pIdx);

    int get_seed_start_pos2(int pIdx);

    int get_seed_end_pos(int pIdx);

    int get_seed_end_pos2(int pIdx);

    int get_length_to_cds(int pIdx);

    void print_all();

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
                           mikan::TRNAStr const &pMiRNA, bool &pEffectiveSite);

    void get_mx_match(mikan::TRNAStr const &pMiRNASeq, mikan::TRNAStr const &pMiRNACompSeq, mikan::TRNAStr const &pMRNASeq,
                      unsigned pSitePos, int pMx, bool &pMatch, bool &pGU, bool &pNoMx, bool &pIsA);

    void get_match_count(mikan::TRNAStr const &pMiRNASeq, mikan::TRNAStr const &pMiRNACompSeq, mikan::TRNAStr const &pMRNASeq,
                         unsigned pSitePos, int pMx1, int pMx2, int &pMatchCount, int &pGUCount);
};


} // namespace tm1p

#endif /* TM1_SEED_SITE_HPP_ */
