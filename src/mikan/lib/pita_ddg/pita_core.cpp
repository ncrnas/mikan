#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"           // MKInput
#include "pita_option.hpp"        // PITAOptions
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_score.hpp"    // PITAMFEScores
#include "pita_site_filter.hpp"   // PITASiteFilter
#include "pita_core.hpp"          // PITACore

namespace ptddg {

//
// PITACore methods
//
int PITACore::write_site_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

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
            score = roundf(score * 100.0f) / 100.0f;

            mOFile1 << toCString(pMiRNAId) << "\t";
            mOFile1 << toCString((seqan::CharString) (mMRNAIds[uniqRNAPosSet[i]])) << "\t";
            mOFile1 << seedStart + 1 << "\t";
            mOFile1 << seedStart + 1 + mikan::SEEDLEN << "\t";
            mOFile1 << toCString((seqan::CharString) (seedTypes[rnaSitePosMap[i][j]])) << "\t";
            mOFile1 << score << "\t";
            mOFile1 << std::endl;

            ++count;
        }

    }

    return 0;

}

int PITACore::write_rna_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &totalScores = mRNAScores.get_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();
    float score;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        score = totalScores[i];
        score = roundf(score * 100.0f) / 100.0f;

        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) (mMRNAIds[mRNAPos[i]])) << "\t";
        mOFile2 << score << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

int PITACore::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    seqan::CharString seedType;
    float dGduplex;
    float dG5;
    float dG3;
    float dGopen;
    float dG0;
    float dG1;
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
            std::cout << "  miRNA:               " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:                "
                      << toCString((seqan::CharString) (mMRNAIds[mRNAPos[uniqRNAPosSet[i]]]));
            std::cout << std::endl;
            std::cout << "  seed type:           " << toCString((seqan::CharString) (seedTypes[rnaSitePosMap[i][j]]))
                      << std::endl;
            std::cout << "  position(start):     " << seedStart + 1 << std::endl;
            std::cout << "  position(end):       " << seedStart + 1 + mikan::SEEDLEN << std::endl;
            std::cout << "  ddG:                 " << score << std::endl;
            std::cout << "  dGduplex(dG5 + dG3): " << dGduplex << std::endl;
            std::cout << "  dG5:                 " << dG5 << std::endl;
            std::cout << "  dG3:                 " << dG3 << std::endl;
            std::cout << "  dGopen(dG0 - dG1):   " << dGopen << std::endl;
            std::cout << "  dG0:                 " << dG0 << std::endl;
            std::cout << "  dG1:                 " << dG1 << std::endl;
            std::cout << std::endl;
        }

    }

    return 0;

}

} // namespace ptddg

