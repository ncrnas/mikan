#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                 // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <map>                   // multimap
#include <utility>               // pair
#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"          // MKInput
#include "ts5_seed_site.hpp"     // TS5SeedSites
#include "ts5_feature.hpp"       // TS5RawFeatures
#include "ts5_site_score.hpp"    // TS5SiteScores
#include "ts5_core.hpp"          // TS5Core

namespace ts5cs {

//
// TS5Core methods
//
void TS5Core::write_site_score_tab(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    int seedStart = sitePos[pSitePosIdx];
    if (seedTypes[pSitePosIdx] == "7mer-A1") {
        seedStart += 1;
    }

    int seedEnd = seedStart + 6;
    if (seedTypes[pSitePosIdx] == "8mer") {
        seedEnd += 1;
    }

    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart << "\t";
    mOFile1 << seedEnd << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << mSiteScores.get_score(pSitePosIdx) << "\t";
    mOFile1 << std::endl;

}

void TS5Core::write_site_score_gff(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    int seedStart = sitePos[pSitePosIdx];
    if (seedTypes[pSitePosIdx] == "7mer-A1") {
        seedStart += 1;
    }

    int seedEnd = seedStart + 6;
    if (seedTypes[pSitePosIdx] == "8mer") {
        seedEnd += 1;
    }

    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart << "\t";
    mOFile1 << seedEnd << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << mSiteScores.get_score(pSitePosIdx);
    mOFile1 << std::endl;

}

void TS5Core::write_rna_score_tab(mikan::TCharStr const &pMiRNAId) {

    const seqan::String<float> &totalScores = mRNAScores.get_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();
    typedef std::multimap<double, unsigned>::iterator TItMap;
    typedef std::pair<double, unsigned> TPosPair;
    TItMap itPos;
    std::multimap<double, unsigned> sortedMRNAByScore;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        sortedMRNAByScore.insert(TPosPair((float) totalScores[i], i));
    }

    for (itPos = sortedMRNAByScore.begin(); itPos != sortedMRNAByScore.end(); ++itPos) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[(*itPos).second]]) << "\t";
        mOFile2 << siteNum[(*itPos).second] << "\t";
        mOFile2 << totalScores[(*itPos).second] << "\t";
        mOFile2 << std::endl;
    }
}

int TS5Core::write_alignment(mikan::TCharStr const &pMiRNAId) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const TS5Alignment &alignment = mSiteScores.get_alignment();
    mikan::TCharStr seedType;
    int seedStart, seedEnd;
    int count = 0;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mSiteScores.mEffectiveSites[i]) {
            continue;
        }

        seedStart = sitePos[i];
        if (seedTypes[i] == "7mer-A1") {
            seedStart += 1;
        }

        seedEnd = seedStart + 6;
        if (seedTypes[i] == "8mer") {
            seedEnd += 1;
        }

        std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
        alignment.write_alignment(i);
        std::cout << "  miRNA:                " << toCString(pMiRNAId) << std::endl;
        std::cout << "  mRNA:                 " << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[i]]) << std::endl;
        std::cout << "  seed type:            " << seedTypes[i] << std::endl;
        std::cout << "  position(seed start): " << seedStart << std::endl;
        std::cout << "  position(seed end):   " << seedEnd << std::endl;
        std::cout << "  context score:        " << mSiteScores.get_score(i);
        std::cout << std::endl << std::endl;

        ++count;

    }

    return 0;
}

} // namespace ts5cs
