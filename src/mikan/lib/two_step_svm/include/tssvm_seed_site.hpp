#ifndef TSSVM_SEED_SITE_HPP_
#define TSSVM_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"          // MKSeedSeqs
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_site_filter.hpp"       // MKSiteFilter
#include "mk_option.hpp"            // MKOptions

namespace tssvm {

//
// Generate miRNA seeds
//
class TSSVMSeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    TSSVMSeedSeqs(mikan::MKOptions const &opts) : MKSeedSeqs(opts) {
        set_flags();
    }

    // Method prototype
    void set_flags();
};

//
// miRNA seed sites
//
class TSSVMSeedSites : public mikan::MKSeedSites {
public:
    // Define methods
    TSSVMSeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {

        mMinToCDS = 15;
        mMinToEnd = 0;

        mUpdatePos = false;
    }

private:
    // Method prototypes
    virtual bool check_position_1(unsigned pMRNAPos, unsigned pSitePos, seqan::CharString &pSeedType);

    virtual bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                           mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                           seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

};

} // namespace tssvm

#endif /* TSSVM_SEED_SITE_HPP_ */
