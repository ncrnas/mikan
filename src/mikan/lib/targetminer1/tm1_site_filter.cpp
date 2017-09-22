#include <seqan/seq_io.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_seed_site.hpp"      // TM1SeedSites
#include "tm1_site_filter.hpp"   // TM1SiteFilter

using namespace seqan;

namespace tm1p {
//
// TM1SiteFilter methods
//
float TM1SiteFilter::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &) {

    mikan::TCharSet const &seedTypes = pSeedSites.get_seed_types();
    mikan::TCharStr seedType = seedTypes[pSitePos];
    float preced;

    if (seedType == "8mer") {
        preced = 0;
    } else if (seedType == "7mer-m8") {
        preced = 1;
    } else if (seedType == "7mer-A1") {
        preced = 2;
    } else if (seedType == "6mer") {
        preced = 3;
    } else {
        preced = 4;
    }

    return preced;
}

} // namespace tm1p
