#ifndef RH2_SEED_SITE_HPP_
#define RH2_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"        // MKSeedSeqs
#include "mk_seed_site.hpp"       // MKSeedSites
#include "mk_option.hpp"          // MKOptions

namespace rh2mfe {

//
// Generate miRNA seeds
//
class RH2SeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    RH2SeedSeqs(mikan::MKOptions const &opts) : MKSeedSeqs(opts) {
        init_from_args();
        set_flags();
    }

    // Method prototypes
    void init_from_args();

    void set_flags();
};

//
// miRNA seed sites
//
class RH2SeedSites : public mikan::MKSeedSites {
public:
    // Define methods
    explicit RH2SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {

        mMinToCDS = 1;
    }

private:
    virtual bool set_new_seed_type(unsigned pMRNAPos, unsigned pSitePos,
                                   mikan::TRNAStr &pMiRNASeq, mikan::TCharSet &pSeedTypeDef,
                                   seqan::CharString &pSeedType, int pMisMatchPos, bool pEffectiveSite);

};

} // namespace rh2mfe

#endif /* RH2_SEED_SITE_HPP_ */
