#ifndef RH2_SEED_SITE_HPP_
#define RH2_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>

namespace rh2mfe{

//
// Store miRNA and mRNA sequences and ids
//
template <class TRNAString>
class RH2Sequences
{
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;

public:
    // Define methods
    RH2Sequences(): mMaxLen(0) {}
    unsigned get_length() const {return length(mSeqIds);}
    seqan::CharString const& get_seq_id(unsigned const pIdx) const {return mSeqIds[pIdx];}
    TRNAString const& get_seq(unsigned const pIdx) const {return mSeqs[pIdx];}
    TCharSet const& get_ids() const {return mSeqIds;}
    TRNASet const& get_seqs() const {return mSeqs;}
    int get_max_seq_len() {return mMaxLen;}

    // Method prototypes
    int read_fasta(seqan::CharString const &pFasta);

private:
    TCharSet mSeqIds;
    TRNASet mSeqs;
    int mMaxLen;
};

//
// Generate miRNA seeds
//
template <class TRNAString>
class RH2SeedSeqs
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Define variables
    seqan::String<bool> mEffectiveSeeds;

public:
    // Define methods
    RH2SeedSeqs() {}
    TRNAString const& get_seed_seq(int i) const {return mSeedSeqs[i];}
    seqan::CharString const& get_seed_type(int i) const {return mSeedTypes[i];}

    // Method prototypes
    int create_seed_seqs(seqan::CharString &pSeedType, seqan::CharString &pOverlapDef);
    void set_mirna_seq(TRNAString pSeq);

private:
    TRNASet mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    TRNAString mMiRNASeq;

private:
    int create_nmer_seed_seqs(TRNAString &pSeedSeq, seqan::CharString &pSeedDef);
    int create_single_guwobble_seed_seqs(TRNAString& pSeedSeq, seqan::CharString &pGUT, seqan::CharString &pGUM);
    int create_multi_guwobble_seed_seqs(TRNAString& pSeedSeq, seqan::CharString &pGUT, seqan::CharString &pGUM);
    int check_redundant_seeds(seqan::CharString &pOverlapDef);

};

//
// miRNA seed sites
//
template <class TRNAString>
class RH2SeedSites
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
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
    RH2SeedSites(TIndexQGram& pRNAIdx, TFinder& pFinder, TRNASet const &pMRNASeqs):
        mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}
    unsigned get_length() const {return seqan::length(mSitePos);}
    seqan::String<unsigned> const& get_mrna_pos() const {return mMRNAPos;}
    seqan::String<unsigned> const& get_site_pos() const {return mSitePos;}
    seqan::StringSet<seqan::CharString> const& get_seed_types() const {return mSeedTypes;}

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
    TIndexQGram& mRNAIdx;
    TFinder& mFinder;
    TRNASet const &mMRNASeqs;

};

} // namespace rh2mfe

#endif /* RH2_SEED_SITE_HPP_ */
