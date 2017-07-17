#ifndef TSSVM_SEED_SITE_HPP_
#define TSSVM_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"          // MKSeedSeqs
#include "mk_seed_site.hpp"         // MKSeedSites

namespace tssvm {

//
// Generate miRNA seeds
//
class TSSVMSeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    TSSVMSeedSeqs() : MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet &pSeedTypeDef);
};

//
// miRNA seed sites
//
class TSSVMSeedSites : public mikan::MKSeedSites {
public:
    // Constant values
    static const unsigned MIN_DIST_TO_CDS = 15;
    static const unsigned MIN_DIST_UTR_END = 0;

    // Define methods
    TSSVMSeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {}

    seqan::String<unsigned> const &get_site_pos_s1() const { return mS1Pos; }

    seqan::String<unsigned> const &get_site_pos_s8() const { return mS8Pos; }

    virtual void clear_pos();

private:
    seqan::String<unsigned> mS1Pos;
    seqan::String<unsigned> mS8Pos;

private:
    virtual bool check_position(unsigned pMRNAPos, unsigned pSitePos, seqan::CharString &pSeedType);

    virtual bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                                   mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                                   seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

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
