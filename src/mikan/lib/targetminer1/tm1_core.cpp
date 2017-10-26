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
void TM1Core::prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    std::string miRNAName = toCString(pMiRNAId);
    std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]);
    unsigned startPos = sitePos[pSitePosIdx] + 1;
    unsigned endPos = sitePos[pSitePosIdx] + 7;
    std::string seedType = toCString((mikan::TCharStr) seedTypes[pSitePosIdx]);
    std::string score1Name = "not used";
    std::string score1 = "0";
    std::string score2Name  = "not used";
    std::string score2 = "0";

    if (mOpts.mGff) {
    } else {
        write_site_score_tab(miRNAName, mRNAName, startPos, endPos, seedType, score1Name, score1, score2Name, score2);
    }

}

void TM1Core::prepare_rna_output(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<float> &scores = mRNAScores.get_scores();
    const seqan::String<int> &predictions = mRNAScores.get_labels();
    const seqan::String<unsigned> &siteNum = mRNAScores.get_site_num();
    mikan::TMRNAPosSet const &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();
    typedef std::multimap<double, unsigned>::iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;
    TItMap itPos;
    std::multimap<double, unsigned> sortedMRNAByScore;

    std::string miRNAName = toCString(pMiRNAId);
    std::string score1Name = "discriminant value";
    std::string score2Name = "predicted label";

    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); ++i) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }
        sortedMRNAByScore.insert(TPosPair(-1.0f * scores[i], i));
    }

    for (itPos = sortedMRNAByScore.begin(); itPos != sortedMRNAByScore.end(); ++itPos) {
        std::stringstream s1, s2;
        std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[uniqRNAPosSet[(*itPos).second]]);
        s1 << scores[(*itPos).second];
        std::string score1 = s1.str();
        s2 << predictions[(*itPos).second];
        std::string score2 = s2.str();

        if (mOpts.mGff) {
        } else {
            write_rna_score_tab(miRNAName, mRNAName, siteNum[(*itPos).second], score1Name, score1, score2Name, score2);
        }

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
