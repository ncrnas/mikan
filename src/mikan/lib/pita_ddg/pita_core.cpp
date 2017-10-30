#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"           // MKInput
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_score.hpp"    // PITAMFEScores
#include "pita_site_filter.hpp"   // PITASiteFilter
#include "pita_core.hpp"          // PITACore

namespace ptddg {

//
// PITACore methods
//
void PITACore::prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    std::stringstream s1;

    std::string miRNAName = toCString(pMiRNAId);
    std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]);
    unsigned startPos = sitePos[pSitePosIdx] + 1;
    unsigned endPos = sitePos[pSitePosIdx] + 7;
    std::string seedType = toCString((mikan::TCharStr) seedTypes[pSitePosIdx]);
    std::string score1Name = "ddG";
    s1 << roundf(mSiteScores.get_score(pSitePosIdx) * 100.0f) / 100.0f;
    std::string score1 = s1.str();
    std::string score2Name  = "not used";
    std::string score2 = "0";

    if (mOpts.mGff) {
        std::string header = "mk-pita site positions and scores - score: alignment, score2: not used";
        std::string src = "mk-pita";
        write_site_score_gff(header, src, miRNAName, mRNAName, startPos, endPos, seedType, score1, score2);
    } else {
        write_site_score_tab(miRNAName, mRNAName, startPos, endPos, seedType, score1Name, score1, score2Name, score2);
    }

}

void PITACore::prepare_rna_output(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<float> &totalScores = mRNAScores.get_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();

    std::string miRNAName = toCString(pMiRNAId);
    std::string score1Name = "ddG";
    std::string score2Name = "not used";

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        std::stringstream s1;
        std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[mRNAPos[i]]);
        s1 << roundf(totalScores[i] * 100.0f) / 100.0f;
        std::string score1 = s1.str();
        std::string score2 = "0";

        if (mOpts.mGff) {
            std::string header = "mk-pita mRNA level scores - score: ddG, score2: not used";
            std::string src = "mk-pita";
            unsigned seq_len = seqan::length(mMRNASeqs[mRNAPos[i]]);
            write_rna_score_gff(header, src, miRNAName, mRNAName, seq_len, siteNum[i], score1, score2);
        } else {
            write_rna_score_tab(miRNAName, mRNAName, siteNum[i], score1Name, score1, score2Name, score2);
        }

    }
}

int PITACore::write_alignment(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    unsigned padw = 22;
    mikan::TCharStr seedType;
    float dGduplex;
    float dG5;
    float dG3;
    float dGopen;
    float dG0;
    float dG1;
    int count = 0;
    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); i++) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        int seedStart;
        float score;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!mSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            seedStart = sitePos[rnaSitePosMap[i][j]];
            score = mSiteScores.get_score(rnaSitePosMap[i][j]);
            score = roundf(score * 100.0f) / 100.0f;

            dGduplex = (float) mSiteScores.get_dgall(rnaSitePosMap[i][j]);
            dGduplex = roundf(dGduplex * 100.0f) / 100.0f;
            dG5 = (float) mSiteScores.get_dg5(rnaSitePosMap[i][j]);
            dG5 = roundf(dG5 * 100.0f) / 100.0f;
            dG3 = (float) mSiteScores.get_dg3(rnaSitePosMap[i][j]);
            dG3 = roundf(dG3 * 100.0f) / 100.0f;
            dGopen =
                    (float) mSiteScores.get_dg0(rnaSitePosMap[i][j]) - (float) mSiteScores.get_dg1(rnaSitePosMap[i][j]);
            dGopen = roundf(dGopen * 100.0f) / 100.0f;
            dG0 = (float) mSiteScores.get_dg0(rnaSitePosMap[i][j]);
            dG0 = roundf(dG0 * 100.0f) / 100.0f;
            dG1 = (float) mSiteScores.get_dg1(rnaSitePosMap[i][j]);
            dG1 = roundf(dG1 * 100.0f) / 100.0f;

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.print_alignment(rnaSitePosMap[i][j]);
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "miRNA: " << toCString(pMiRNAId) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "mRNA: " << toCString((mikan::TCharStr) (mMRNAIds[mRNAPos[uniqRNAPosSet[i]]])) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "seed type: " << toCString((mikan::TCharStr) (seedTypes[rnaSitePosMap[i][j]])) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "start (1-base): " << seedStart + 1 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "end (1-base): " << seedStart + 1 + mikan::SEEDLEN << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "ddG: " << score << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "dGduplex(dG5 + dG3): " << dGduplex << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "dG5: " << dG5 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "dG3: " << dG3 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "dGopen(dG0 - dG1): " << dGopen << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "dG0: " << dG0 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "dG1: " << dG1 << std::endl;
            std::cout << std::endl;

            ++count;
        }

    }

    return 0;

}

} // namespace ptddg

