#ifndef TSSVM_SEED_SITE_HPP_
#define TSSVM_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder

namespace tssvm {

//
// Generate miRNA seeds
//
class TSSVMSeedSeqs {
public:
    // Define variables
    seqan::String<bool> mEffectiveSeeds;

public:
    // Define methods
    TSSVMSeedSeqs() {
        resize(mRNAChar, 4);
        mRNAChar[0] = 'A';
        mRNAChar[1] = 'C';
        mRNAChar[2] = 'G';
        mRNAChar[3] = 'U';
    }

    mikan::TRNAStr const &get_seed_seq(int i) const { return mSeedSeqs[i]; }

    seqan::CharString const &get_seed_type(int i) const { return mSeedTypes[i]; }

    unsigned get_mismatched_pos(int i) { return mMisMatchPos[i]; }

    // Method prototypes
    int create_seed_seqs();

    void set_mirna_seq(mikan::TRNAStr pSeq);

private:
    mikan::TRNASet mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    mikan::TRNAStr mMiRNASeq;
    seqan::String<unsigned> mMisMatchPos;
    mikan::TRNAStr mRNAChar;

private:
    int create_non_stringent_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_guwobble_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_lp_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_bm_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int create_bt_seed_seqs(mikan::TRNAStr &pSeedSeq);

    int add_seeds_in_reverse_order(mikan::TRNASet &pSeedSeqs, seqan::StringSet<seqan::CharString> &pSeedTypes,
                                   seqan::String<unsigned> &pMisMatchPos);

};

//
// miRNA seed sites
//
class TSSVMSeedSites {
public:

    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 15;
    static const unsigned MIN_DIST_UTR_END = 0;

    // Define variables
    seqan::String<bool> mEffectiveSites;

public:
    // Define methods
    TSSVMSeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            mRNAIdx(pRNAIdx), mFinder(pFinder), mMRNASeqs(pMRNASeqs) {}

    unsigned get_length() const { return seqan::length(mSitePos); }

    seqan::String<unsigned> const &get_mrna_pos() const { return mMRNAPos; }

    seqan::String<unsigned> const &get_site_pos() const { return mSitePos; }

    seqan::String<unsigned> const &get_site_pos_s1() const { return mS1Pos; }

    seqan::String<unsigned> const &get_site_pos_s8() const { return mS8Pos; }

    seqan::StringSet<seqan::CharString> const &get_seed_types() const { return mSeedTypes; }

    seqan::String<unsigned> const &get_mismatched_pos() const { return mMisMatchPos; }

    // Method prototypes
    void reset_finder();

    int find_seed_sites(mikan::TRNAStr const &pMiRNA);

    void clear_pos();

private:
    seqan::String<unsigned> mMRNAPos;
    seqan::String<unsigned> mSitePos;
    seqan::StringSet<seqan::CharString> mSeedTypes;
    seqan::String<unsigned> mS1Pos;
    seqan::String<unsigned> mS8Pos;
    seqan::String<unsigned> mMisMatchPos;
    mikan::TIndexQGram &mRNAIdx;
    mikan::TFinder &mFinder;
    mikan::TRNASet const &mMRNASeqs;

private:
    int set_seed_pos(TSSVMSeedSeqs &pSeedSeqs, unsigned pMRNAPos, unsigned pSitePos,
                     mikan::TRNAStr const &pMiRNA, const mikan::TRNAStr &pSeedSeq, unsigned pIdx);

    int set_seed_type(const seqan::CharString &pCurType, const mikan::TRNAStr &pMRNASeq,
                      const mikan::TRNAStr &pMiRNASeq, unsigned pM8Pos, unsigned pA1Pos, unsigned pMisMatchedPos,
                      const mikan::TRNAStr &pSeedSeq);

};

//
// miRNA seed sites overlap
//
class TSSVMSeedSiteOverlap {
public:
    // Define methods
    TSSVMSeedSiteOverlap() {}

    // Method prototype
    int filter_overlapped_sites(TSSVMSeedSites &pSeedSites, unsigned pMRNANum);

    void clear_site_pos();

    std::set<unsigned> &get_mrna_pos_set() { return mRNAPosSet; }

    std::multimap<unsigned, unsigned> &get_site_map() { return mSiteMap; }

    seqan::StringSet<seqan::String<unsigned> > &get_sorted_mrna_pos() { return mSortedMRNAPos; }

private:
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<TItMap, TItMap> TItRetPair;
    typedef std::multimap<unsigned, unsigned>::iterator TITSeedTypes;
    typedef std::multimap<unsigned, unsigned>::iterator TITPos;
    typedef std::pair<unsigned, unsigned> TPosPair;

    std::set<unsigned> mRNAPosSet;
    std::multimap<unsigned, unsigned> mSiteMap;
    seqan::StringSet<seqan::String<unsigned> > mSortedMRNAPos;

private:
    void cluster_site_pos(TSSVMSeedSites &pSeedSites);

    void sort_by_seed_type(TSSVMSeedSites &pSeedSites, int pPosIdx);

    unsigned get_seedtype_precedence(const seqan::CharString &pSeedType);

    void mark_overlapped_sites(TSSVMSeedSites &pSeedSites,
                               std::multimap<unsigned, unsigned> &pSortedSeeds);
};

} // namespace tssvm

#endif /* TSSVM_SEED_SITE_HPP_ */
