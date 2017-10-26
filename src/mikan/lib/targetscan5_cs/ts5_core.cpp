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
void TS5Core::prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    std::stringstream s1;

    std::string miRNAName = toCString(pMiRNAId);
    std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]);
    unsigned startPos = sitePos[pSitePosIdx];
    if (seedTypes[pSitePosIdx] == "7mer-A1") {
        startPos += 1;
    }
    unsigned endPos = startPos + 6;
    if (seedTypes[pSitePosIdx] == "8mer") {
        endPos += 1;
    }
    std::string seedType = toCString((mikan::TCharStr) seedTypes[pSitePosIdx]);
    std::string score1Name = "context score";
    s1 << mSiteScores.get_score(pSitePosIdx);
    std::string score1 = s1.str();
    std::string score2Name  = "not used";
    std::string score2 = "0";

    if (mOpts.mGff) {
    } else {
        write_site_score_tab(miRNAName, mRNAName, startPos, endPos, seedType, score1Name, score1, score2Name, score2);
    }

}

void TS5Core::prepare_rna_output(mikan::TCharStr const &pMiRNAId) {

    const seqan::String<float> &totalScores = mRNAScores.get_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();
    typedef std::multimap<double, unsigned>::iterator TItMap;
    typedef std::pair<double, unsigned> TPosPair;
    TItMap itPos;
    std::multimap<double, unsigned> sortedMRNAByScore;

    std::string miRNAName = toCString(pMiRNAId);
    std::string score1Name = "context score";
    std::string score2Name = "not used";

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        sortedMRNAByScore.insert(TPosPair((float) totalScores[i], i));
    }

    for (itPos = sortedMRNAByScore.begin(); itPos != sortedMRNAByScore.end(); ++itPos) {
        std::stringstream s1;
        std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[mRNAPos[(*itPos).second]]);
        s1 << totalScores[(*itPos).second];
        std::string score1 = s1.str();
        std::string score2 = "0";

        if (mOpts.mGff) {
        } else {
            write_rna_score_tab(miRNAName, mRNAName, siteNum[(*itPos).second], score1Name, score1, score2Name, score2);
        }

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

    unsigned padw = 17;

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
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "miRNA: " << toCString(pMiRNAId) << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "mRNA: " << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[i]]) << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "seed type: " << seedTypes[i] << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "start (1-base): " << seedStart << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "end (1-base): " << seedEnd << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "context score: " << mSiteScores.get_score(i) << std::endl;
        std::cout << std::endl;

        ++count;

    }

    return 0;
}

} // namespace ts5cs
