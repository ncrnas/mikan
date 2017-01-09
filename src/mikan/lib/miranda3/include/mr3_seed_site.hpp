#ifndef MR3_SEED_SITE_HPP_
#define MR3_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>

namespace mr3as{

//
// Store miRNA and mRNA sequences and ids
//
template <class TRNAString>
class MR3Sequences
{
public:
    // Define types
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Define methods
    MR3Sequences(): mMaxLen(0) {}
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
class MR3SeedSeqs
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;

    // Define variables
    seqan::String<bool> mEffectiveSeeds;

public:
    // Define methods
    MR3SeedSeqs() {}
    TRNAString const& get_seed_seq(int i) const {return mSeedSeqs[i];}
    seqan::CharString const& get_seed_type(int i) const {return mSeedTypes[i];}
    unsigned get_mismatched_pos(int i) {return mMisMatchPos[i];}

    // Method prototypes
    int create_seed_seqs(seqan::StringSet<seqan::CharString> &pSeedType);
    void set_mirna_seq(TRNAString pSeq);

private:
    TRNASet mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    seqan::String<unsigned> mMisMatchPos;
    TRNAString mMiRNASeq;

private:
    int create_nmer_seed_seqs(TRNAString &pSeedSeq, seqan::StringSet<seqan::CharString> &pSeedDef);
    int create_single_guwobble_seed_seqs(TRNAString &pSeedSeq);
    int create_multi_guwobble_seed_seqs(TRNAString &pSeedSeq);
    int create_mismatch_seed_seqs(TRNAString &pSeedSeq, bool pIsGUMM=false, int pGUPos=0);
    int create_gu_mismatch_seed_seqs(TRNAString &pSeedSeq);
    int create_bt_seed_seqs(TRNAString &pSeedSeq);
    int check_redundant_seeds();

};

//
// miRNA seed sites
//
template <class TRNAString>
class MR3SeedSites
{
public:
    // Define types
    typedef seqan::StringSet<TRNAString> TRNASet;
    typedef seqan::StringSet<seqan::CharString> TCharSet;
    typedef seqan::Index<TRNASet, seqan::IndexQGram<seqan::UngappedShape<6> > > TIndexQGram;
    typedef seqan::Finder<TIndexQGram> TFinder;

    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 7;
    static const unsigned MIN_DIST_UTR_END = 0;
    static const unsigned INDEXED_SEQ_LEN = 6;
    static const bool FORCE_LAST_MATCH = false;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    MR3SeedSites(TIndexQGram& pRNAIdx, TFinder& pFinder, TRNASet const &pMRNASeqs):
        mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}
    unsigned get_length() const {return seqan::length(mSitePos);}
    seqan::String<unsigned> const& get_mrna_pos() const {return mMRNAPos;}
    seqan::String<unsigned> const& get_site_pos() const {return mSitePos;}
    seqan::StringSet<seqan::CharString> const& get_seed_types() const {return mSeedTypes;}
    seqan::String<int> const& get_mismatched_pos() const {return mMisMatchPos;}

    // Method prototypes
    void reset_finder();
    int find_seed_sites(TRNAString const &pMiRNA, seqan::StringSet<seqan::CharString> &pSeedDef);
    void clear_pos();

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    seqan::String<int> mMisMatchPos;
    TIndexQGram& mRNAIdx;
    TFinder& mFinder;
    TRNASet const &mMRNASeqs;

private:
    void set_new_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
            unsigned pMRNAPos, unsigned pSitePos, TRNAString const &pMiRNA, unsigned pMisMatchPos,
            bool &pEffectiveSite);

    void set_mx_matches(unsigned pMRNAPos, unsigned pSitePos, TRNAString const &pMiRNA, int pMx,
            bool &pMatchMx, bool &pGutMx, bool &pGumMx);

    void set_stringent_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
            bool pMatchMx8, bool pMatchMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);
    void set_single_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
            int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
            bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);
    void set_multiple_gu_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
            int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
            bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);
    void set_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
            int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
            bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);
    void set_gu_mismatch_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
            int pM1, int pM2, bool pMatchMx8, bool pMatchMx9, bool pGutMx8, bool pGutMx9,
            bool pGumMx8, bool pGumMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);
    void set_bt_seed_type(unsigned pMRNAPos, unsigned pSitePos, TRNAString const &pMiRNA, unsigned pMisMatchPos,
            seqan::CharString &pNewSeedType);
    void set_6mer_seed_type(seqan::CharString &pCurSeedType, seqan::StringSet<seqan::CharString> &pSeedDef,
            bool pMatchMx8, bool pMatchMx9, unsigned pMisMatchPos, seqan::CharString &pNewSeedType);

    void check_last_match(bool pMatchM8, bool pMatchM9, seqan::CharString &pNewSeedType);

};

} // namespace mr3as

#endif /* MR3_SEED_SITE_HPP_ */
