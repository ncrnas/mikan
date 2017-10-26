#include <math.h>                // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                 // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"          // MKInput
#include "rh2_seed_site.hpp"     // RH2Sequences, RH2SeedSites
#include "rh2_site_score.hpp"    // RH2SiteScores, RH2TotalScores
#include "rh2_site_filter.hpp"   // RH2Overlap, RH2TopNSites, RH2SortedSitePos
#include "rh2_core.hpp"          // RH2Core

namespace rh2mfe {

//
// RH2Core methods
//
void RH2Core::prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    std::stringstream s1, s2;

    std::string miRNAName = toCString(pMiRNAId);
    std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]);
    unsigned startPos = sitePos[pSitePosIdx] + 1;
    unsigned endPos = sitePos[pSitePosIdx] + 7;
    std::string seedType = toCString((mikan::TCharStr) seedTypes[pSitePosIdx]);
    std::string score1Name = "MFE";
    s1 << roundf(mSiteScores.get_score(pSitePosIdx) * 10.0f) / 10.0f;
    std::string score1 = s1.str();
    std::string score2Name  = "normalized";
    s2 << mSiteScores.get_norm_score(pSitePosIdx);
    std::string score2 = s2.str();

    if (mOpts.mGff) {
    } else {
        write_site_score_tab(miRNAName, mRNAName, startPos, endPos, seedType, score1Name, score1, score2Name, score2);
    }

}

void RH2Core::prepare_rna_output(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<float> &mfeScores = mRNAScores.get_mfe_minlogtotal();
    const seqan::String<float> &normScores = mRNAScores.get_norm_maxlogtotal();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();

    if (mPrintRNAheader) {
        mOFile2 << "# miRNA name, ";
        mOFile2 << "mRNA name, ";
        mOFile2 << "number of sites, ";
        mOFile2 << "score 1 (MFE), ";
        mOFile2 << "score 2 (normalized)";
        mOFile2 << std::endl;
        mPrintRNAheader = false;
    }

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << mfeScores[i] << "\t";
        mOFile2 << normScores[i];
        mOFile2 << std::endl;
    }

}

int RH2Core::write_alignment(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    unsigned padw = 19;

    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); i++) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        int seedStart;
        int count = 0;
        float score;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!mSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            seedStart = sitePos[rnaSitePosMap[i][j]];
            score = mSiteScores.get_score(rnaSitePosMap[i][j]);
            score = roundf(score * 10.0f) / 10.0f;

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.write_alignment(rnaSitePosMap[i][j]);
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "miRNA: " << toCString(pMiRNAId) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "mRNA: " << toCString((mikan::TCharStr) mMRNAIds[uniqRNAPosSet[i]]) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "seed type: " << toCString((mikan::TCharStr) seedTypes[rnaSitePosMap[i][j]])  << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "start (1-base): " << seedStart + 1 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "end (1-base): " << seedStart + 7 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "mfe: " << score << " kcal/mol" << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "normalized score: " << mSiteScores.get_norm_score(rnaSitePosMap[i][j]) << std::endl;
            std::cout << std::endl;

            ++count;
        }

    }

    return 0;

}

} // namespace rh2mfe

