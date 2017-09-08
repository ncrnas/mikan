#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "pita_site_score.hpp"    // PITAMFEScores
#include "pita_site_filter.hpp"   // PITASiteFilter, PITASortedSitePos

using namespace seqan;

namespace ptddg {

//
// PITASiteFilter methods
//
float PITASiteFilter::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &) {

    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    seqan::CharString seedType = seedTypes[pSitePos];
    float preced;

    if (seedType == "8mer") {
        preced = 0;
    } else if (seedType == "7mer") {
        preced = 1;
    } else if (seedType == "6mer") {
        preced = 2;
    } else if (seedType == "8mer_GUT" || seedType == "8mer_GUM") {
        preced = 3;
    } else if (seedType == "8mer_GU+") {
        preced = 4;
    } else if (seedType == "8mer_MM") {
        preced = 5;
    } else if (seedType == "7mer_GUT" || seedType == "7mer_GUM") {
        preced = 6;
    } else if (seedType == "7mer_GU+") {
        preced = 7;
    } else if (seedType == "7mer_MM") {
        preced = 8;
    } else if (seedType == "8mer_MMGU") {
        preced = 9;
    } else if (seedType == "7mer_MMGU") {
        preced = 10;
    } else if (seedType == "8mer_MMGU+") {
        preced = 11;
    } else if (seedType == "7mer_MMGU+") {
        preced = 12;
    } else {
        preced = 13;
    }

    return preced;
}

} // namespace ptddg
