#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "pita_score.hpp"         // PITAMFEScores
#include "pita_site_cluster.hpp"  // PITASiteFilter, PITASortedSitePos

using namespace seqan;

namespace ptddg {

//
// PITASiteFilter methods
//
unsigned PITASiteFilter::get_seedtype_precedence(CharString const &pSeedType) {
    unsigned preced;

    if (pSeedType == "8mer") {
        preced = 0;
    } else if (pSeedType == "7mer") {
        preced = 1;
    } else if (pSeedType == "6mer") {
        preced = 2;
    } else if (pSeedType == "8mer_GUT" || pSeedType == "8mer_GUM") {
        preced = 3;
    } else if (pSeedType == "8mer_GU+") {
        preced = 4;
    } else if (pSeedType == "8mer_MM") {
        preced = 5;
    } else if (pSeedType == "7mer_GUT" || pSeedType == "7mer_GUM") {
        preced = 6;
    } else if (pSeedType == "7mer_GU+") {
        preced = 7;
    } else if (pSeedType == "7mer_MM") {
        preced = 8;
    } else if (pSeedType == "8mer_MMGU") {
        preced = 9;
    } else if (pSeedType == "7mer_MMGU") {
        preced = 10;
    } else if (pSeedType == "8mer_MMGU+") {
        preced = 11;
    } else if (pSeedType == "7mer_MMGU+") {
        preced = 12;
    } else {
        preced = 13;
    }

    return preced;
}

void PITASiteFilter::set_intervals(
        mikan::MKSeedSites &pSeedSites,
        unsigned pSiteIdx,
        unsigned &pStartSearch,
        unsigned &pEndSearch,
        unsigned &pStartAdd,
        unsigned &pEndAdd,
        bool &pSearchOverlap) {

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();

    pStartSearch = sitePos[pSiteIdx];
    pEndSearch = pStartSearch + mGapLen;
    pStartAdd = pStartSearch;
    pEndAdd = pEndSearch;
    pSearchOverlap = true;
}

} // namespace ptddg
