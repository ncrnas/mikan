#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"           // MKInput
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "mr3_site_score.hpp"     // MR3SiteScores
#include "mr3_site_filter.hpp"    // MR3SiteFilter
#include "mr3_core.hpp"           // MR3Core

namespace mr3as {

//
// MR3Core methods
//
void MR3Core::prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    std::stringstream s1, s2;

    std::string miRNAName = toCString(pMiRNAId);
    std::string mRNAName = toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]);
    unsigned startPos = sitePos[pSitePosIdx] + 1;
    unsigned endPos = sitePos[pSitePosIdx] + 7;
    std::string seedType = toCString((mikan::TCharStr) seedTypes[pSitePosIdx]);
    std::string score1Name = "alignment";
    s1 << roundf(mSiteScores.get_align_score(pSitePosIdx) * 100.0f) / 100.0f;
    std::string score1 = s1.str();
    std::string score2Name  = "MFE";
    s2 << roundf(mSiteScores.get_energy_score(pSitePosIdx) * 100.0f) / 100.0f;
    std::string score2 = s2.str();

    if (mOpts.mGff) {
    } else {
        write_site_score_tab(miRNAName, mRNAName, startPos, endPos, seedType, score1Name, score1, score2Name, score2);
    }

}

void MR3Core::prepare_rna_output(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<float> &alignScores = mRNAScores.get_align_maxlogtotal();
    const seqan::String<float> &enScores = mRNAScores.get_energy_minlogtotal();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();

    if (mPrintRNAheader) {
        mOFile2 << "# miRNA name, ";
        mOFile2 << "mRNA name, ";
        mOFile2 << "number of sites, ";
        mOFile2 << "score 1 (alignment), ";
        mOFile2 << "score 2 (MFE)";
        mOFile2 << std::endl;
        mPrintRNAheader = false;
    }

    for (unsigned i = 0; i < length(mRNAPos); ++i) {

        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << alignScores[i] << "\t";
        mOFile2 << enScores[i];
        mOFile2 << std::endl;
    }
}

int MR3Core::write_alignment(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    unsigned padw = 18;
    unsigned count = 0;
    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); i++) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!mSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            mikan::TCharStr seedType = seedTypes[rnaSitePosMap[i][j]];
            int seedStart = sitePos[rnaSitePosMap[i][j]];
            float align_score = mSiteScores.get_align_score(rnaSitePosMap[i][j]);
            align_score = roundf(align_score * 100.0f) / 100.0f;
            float energy_score = mSiteScores.get_energy_score(rnaSitePosMap[i][j]);
            energy_score = roundf(energy_score * 100.0f) / 100.0f;

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.print_alignment(rnaSitePosMap[i][j]);
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "miRNA: " << toCString(pMiRNAId) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "mRNA: " << toCString((mikan::TCharStr) mMRNAIds[uniqRNAPosSet[i]]) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "seed type: " << toCString((mikan::TCharStr) seedTypes[rnaSitePosMap[i][j]]) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "start (1-base): " << seedStart + 1 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "end (1-base): " << seedStart + 1 + mikan::SEEDLEN << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "alignment score: " << align_score << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "energy score: " << energy_score << std::endl;
            std::cout << std::endl;

            ++count;
        }
    }

    return 0;

}

} // namespace mr3as
