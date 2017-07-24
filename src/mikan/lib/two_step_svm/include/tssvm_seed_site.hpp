#ifndef TSSVM_SEED_SITE_HPP_
#define TSSVM_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"          // MKSeedSeqs
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_site_filter.hpp"       // MKSiteFilter

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
    // Define methods
    TSSVMSeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {
        mCheckPosMethod = "tssvm";
        mMinToCDS = 15;
        mMinToEnd = 0;

        mUpdatePos = false;
    }

private:
    bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                           mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                           seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

};

//
// miRNA seed sites overlap
//
class TSSVMSiteFilter : public mikan::MKSiteFilter {
public:
    // Define types
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<unsigned, unsigned> TPosPair;

    // Define methods
    TSSVMSiteFilter() : MKSiteFilter() {}

private:
    unsigned get_seedtype_precedence(seqan::CharString const &pSeedType);

    void set_intervals(mikan::MKSeedSites &pSeedSites, unsigned pSiteIdx, unsigned &pStartSearch,
                       unsigned &pEndSearch, unsigned &pStartAdd, unsigned &pEndAdd,
                       bool &pSearchOverlap);

};

} // namespace tssvm

#endif /* TSSVM_SEED_SITE_HPP_ */
