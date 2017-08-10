#ifndef MKE_SEED_SITE_HPP_
#define MKE_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"       // MKSeedSeqs
#include "mk_seed_site.hpp"      // MKSeedSites
#include "mk_option.hpp"         // MKOptions

namespace mkens {

//
// Generate miRNA seeds
//
class MKESeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    MKESeedSeqs(mikan::MKOptions const &opts) : MKSeedSeqs(opts) {
        set_flags();
    }

    // Method prototype
    void set_flags();
};

//
// miRNA seed sites
//
class MKESeedSites : public mikan::MKSeedSites {
public:
    // Define methods
    MKESeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {}

};

} // namespace mkens

#endif /* MKE_SEED_SITE_HPP_ */
