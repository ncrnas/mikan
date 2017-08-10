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
int MR3Core::write_site_score(seqan::CharString const &pMiRNAId) {

    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); i++) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        int seedStart;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!mSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            seedStart = sitePos[rnaSitePosMap[i][j]];
            float scoreAlign = mSiteScores.get_align_score(rnaSitePosMap[i][j]);
            scoreAlign = roundf(scoreAlign * 100.0f) / 100.0f;
            float scoreEn = mSiteScores.get_energy_score(rnaSitePosMap[i][j]);
            scoreEn = roundf(scoreEn * 100.0f) / 100.0f;

            mOFile1 << toCString(pMiRNAId) << "\t";
            mOFile1 << toCString((seqan::CharString) mMRNAIds[uniqRNAPosSet[i]]) << "\t";
            mOFile1 << seedStart + 1 << "\t";
            mOFile1 << seedStart + 1 + mikan::SEEDLEN << "\t";
            mOFile1 << toCString((seqan::CharString) seedTypes[rnaSitePosMap[i][j]]) << "\t";
            mOFile1 << scoreAlign << "\t";
            mOFile1 << scoreEn << "\t";
            mOFile1 << std::endl;
        }

    }

    return 0;

}

int MR3Core::write_rna_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &totalAlignScores = mRNAScores.get_align_scores();
    const seqan::String<float> &totalEnScores = mRNAScores.get_energy_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {

        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile2 << totalAlignScores[i] << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << totalEnScores[i] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

int MR3Core::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); i++) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!mSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            seqan::CharString seedType = seedTypes[rnaSitePosMap[i][j]];
            int seedStart = sitePos[rnaSitePosMap[i][j]];
            float align_score = mSiteScores.get_align_score(rnaSitePosMap[i][j]);
            align_score = roundf(align_score * 100.0f) / 100.0f;
            float energy_score = mSiteScores.get_energy_score(rnaSitePosMap[i][j]);
            energy_score = roundf(energy_score * 100.0f) / 100.0f;

            std::cout << "### " << (i + j) + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.print_alignment(rnaSitePosMap[i][j]);
            std::cout << "  miRNA:           " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:            " << toCString((seqan::CharString) mMRNAIds[uniqRNAPosSet[i]])
                      << std::endl;
            std::cout << "  seed type:       " << toCString((seqan::CharString) seedTypes[rnaSitePosMap[i][j]])
                      << std::endl;
            std::cout << "  position(start): " << seedStart + 1 << std::endl;
            std::cout << "  position(end):   " << seedStart + 1 + mikan::SEEDLEN << std::endl;
            std::cout << "  alignment score: " << align_score << std::endl;
            std::cout << "  energy score:    " << energy_score << std::endl;

            std::cout << std::endl;


        }
    }

    return 0;

}

} // namespace mr3as
