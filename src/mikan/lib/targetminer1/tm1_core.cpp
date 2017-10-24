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
#include "tm1_core.hpp"          // TM1Core

namespace tm1p {

//
// TM1Core methods
//
void TM1Core::write_site_score_tab(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    int seedStart = sitePos[pSitePosIdx] + 1;
    int seedEnd = seedStart + 6;

    if (mPrintSiteHeader) {
        mOFile1 << "# miRNA name, ";
        mOFile1 << "mRNA name, ";
        mOFile1 << "start (1-base), ";
        mOFile1 << "end (1-base), ";
        mOFile1 << "seed type, ";
        mOFile1 << "score 1 (not used), ";
        mOFile1 << "score 2 (not used)";
        mOFile1 << std::endl;
        mPrintSiteHeader = false;
    }

    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart << "\t";
    mOFile1 << seedEnd << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << 0 << "\t";
    mOFile1 << 0;
    mOFile1 << std::endl;

}

void TM1Core::write_site_score_gff(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    int seedStart = sitePos[pSitePosIdx];
    int seedEnd = seedStart + 6;

    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart << "\t";
    mOFile1 << seedEnd << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << 0;
    mOFile1 << std::endl;

}

void TM1Core::write_rna_score_tab(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<float> &scores = mRNAScores.get_scores();
    const seqan::String<int> &predictions = mRNAScores.get_labels();
    const seqan::String<unsigned> &siteNum = mRNAScores.get_site_num();
    mikan::TMRNAPosSet const &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();
    typedef std::multimap<double, unsigned>::iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;
    TItMap itPos;
    std::multimap<double, unsigned> sortedMRNAByScore;

    if (mPrintRNAheader) {
        mOFile2 << "# miRNA name, ";
        mOFile2 << "mRNA name, ";
        mOFile2 << "number of sites, ";
        mOFile2 << "score 1 (discriminant value), ";
        mOFile2 << "score 2 (predicted label)";
        mOFile2 << std::endl;
        mPrintRNAheader = false;
    }

    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); ++i) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }
        sortedMRNAByScore.insert(TPosPair(-1.0f * scores[i], i));
    }

    for (itPos = sortedMRNAByScore.begin(); itPos != sortedMRNAByScore.end(); ++itPos) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((mikan::TCharStr) mMRNAIds[uniqRNAPosSet[(*itPos).second]]) << "\t";
        mOFile2 << siteNum[(*itPos).second] << "\t";
        mOFile2 << scores[(*itPos).second] << "\t";
        mOFile2 << predictions[(*itPos).second];
        mOFile2 << std::endl;
    }
}

int TM1Core::write_alignment(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &mSeedTypes = mSeedSites.get_seed_types();

    unsigned padw = 17;

    mikan::TCharStr seedType;
    int seedStart, seedEnd;
    int count = 0;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mSeedSites.mEffectiveSites[i]) {
            continue;
        }

        seedType = mSeedTypes[i];
        seedStart = sitePos[i] + 1;
        seedEnd = seedStart + 6;

        std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
        mSiteScores.write_alignment(i);
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "miRNA: " << toCString(pMiRNAId) << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "mRNA: " << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[i]]) << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "seed type: " << toCString(seedType) << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "start (1-base): " << seedStart << std::endl;
        std::cout << std::right << std::setw(padw) << std::setfill(' ');
        std::cout << "end (1-base): " << seedEnd << std::endl;
        std::cout << std::endl;

        ++count;

    }

    return 0;
}

} // namespace tm1p
