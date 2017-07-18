#ifndef TS5_SEED_SITE_HPP_
#define TS5_SEED_SITE_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"       // MKSeedSeqs
#include "mk_seed_site.hpp"      // MKSeedSites

namespace ts5cs {

//
// Generate miRNA seeds
//
class TS5SeedSeqs : public mikan::MKSeedSeqs {
public:
    // Define methods
    TS5SeedSeqs(): MKSeedSeqs() {}

    // Method prototypes
    void set_flags(mikan::TCharSet &pSeedTypeDef);
};

//
// miRNA seed sites
//
class TS5SeedSites : public mikan::MKSeedSites {
public:
    // Define methods
    TS5SeedSites(mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder, mikan::TRNASet const &pMRNASeqs) :
            MKSeedSites(pRNAIdx, pFinder, pMRNASeqs) {}

};

} // namespace ts5cs

#endif /* TS5_SEED_SITE_HPP_ */
