#include <seqan/seq_io.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_site_filter.hpp"    // TSSVMSiteFilter

using namespace seqan;

namespace tssvm {

//
// TSSVMSiteFilter methods
//
float TSSVMSiteFilter::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &) {

    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();
    seqan::CharString seedType = seedTypes[pSitePos];
    float preced;

    if (seedType == "8mer" || seedType == "7mer-A1" || seedType == "7mer-m8") {
        preced = 0;
    } else if (seedType == "6mer") {
        preced = 1;
    } else if (seedType == "GUM") {
        preced = 2;
    } else if (seedType == "GUT") {
        preced = 3;
    } else if (seedType == "BT") {
        preced = 4;
    } else if (seedType == "BM") {
        preced = 5;
    } else if (seedType == "LP") {
        preced = 6;
    } else {
        preced = 7;
    }

    return preced;
}

void TSSVMSiteFilter::set_intervals(
        mikan::MKSeedSites &pSeedSites,
        mikan::MKSiteScores &,
        unsigned pSiteIdx,
        unsigned &pStartSearch,
        unsigned &pEndSearch,
        unsigned &pStartAdd,
        unsigned &pEndAdd,
        bool &pSearchOverlap) {

    String<unsigned> const &sitePos = pSeedSites.get_site_pos();
    StringSet<CharString> const &seedTypes = pSeedSites.get_seed_types();

    pStartSearch = 0;
    pEndSearch = 0;
    pStartAdd = sitePos[pSiteIdx] - 1;
    pEndAdd = sitePos[pSiteIdx] + 7;
    pSearchOverlap = true;

    if (seedTypes[pSiteIdx] == "8mer"
        || seedTypes[pSiteIdx] == "7mer-A1"
        || seedTypes[pSiteIdx] == "7mer-m8"
        || seedTypes[pSiteIdx] == "6mer") {
        pStartSearch = pStartAdd;
        pEndSearch = pEndAdd;
        pSearchOverlap = false;
    } else if (seedTypes[pSiteIdx] == "GUM"
               || seedTypes[pSiteIdx] == "GUT"
               || seedTypes[pSiteIdx] == "LP") {
        pStartSearch = pStartAdd - 1;
        pEndSearch = pEndAdd + 1;
    } else if (seedTypes[pSiteIdx] == "BT") {
        pStartSearch = pStartAdd - 1;
        pEndSearch = pEndAdd + 2;
    } else if (seedTypes[pSiteIdx] == "BM") {
        pStartAdd = sitePos[pSiteIdx];
        pEndAdd = sitePos[pSiteIdx] + 8;
        pStartSearch = pStartAdd - 1;
        pEndSearch = pEndAdd;
    }
}

} // namespace tssvm
