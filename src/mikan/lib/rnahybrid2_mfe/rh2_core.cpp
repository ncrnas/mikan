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
void RH2Core::write_site_score_tab(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    int seedStart = sitePos[pSitePosIdx];
    float score = mSiteScores.get_score(pSitePosIdx);
    score = roundf(score * 10.0f) / 10.0f;

    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart + 1 << "\t";
    mOFile1 << seedStart + 7 << "\t";
    //        mOFile1 << mSiteScores.get_wide_site_start(posIdx) + 1  << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << score << "\t";
    mOFile1 << mSiteScores.get_norm_score(pSitePosIdx);
    mOFile1 << std::endl;

}

void RH2Core::write_site_score_gff(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    int seedStart = sitePos[pSitePosIdx];
    float score = mSiteScores.get_score(pSitePosIdx);
    score = roundf(score * 10.0f) / 10.0f;

    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart + 1 << "\t";
    mOFile1 << seedStart + 7 << "\t";
    //        mOFile1 << mSiteScores.get_wide_site_start(posIdx) + 1  << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << score << "\t";
    mOFile1 << mSiteScores.get_norm_score(pSitePosIdx);
    mOFile1 << std::endl;

}

void RH2Core::write_rna_score_tab(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<float> &mfeScores = mRNAScores.get_mfe_minlogtotal();
    const seqan::String<float> &normScores = mRNAScores.get_norm_maxlogtotal();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile2 << mfeScores[i] << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << normScores[i] << "\t";
        mOFile2 << std::endl;
    }

}

int RH2Core::write_alignment(mikan::TCharStr const &pMiRNAId) {
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();

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
            score = roundf(score * 10.0f) / 10.0f;

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.write_alignment(rnaSitePosMap[i][j]);
            std::cout << "  miRNA:               " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:                " << toCString((mikan::TCharStr) mMRNAIds[uniqRNAPosSet[i]])
                      << std::endl;
            std::cout << "  seed type:           " << toCString((mikan::TCharStr) seedTypes[rnaSitePosMap[i][j]])
                      << std::endl;
            std::cout << "  position(target 5'): " << mSiteScores.get_wide_site_start(rnaSitePosMap[i][j]) + 1;
            std::cout << std::endl;
            std::cout << "  position(seed):      " << seedStart + 1 << std::endl;
            std::cout << "  mfe:                 " << score << " kcal/mol" << std::endl;
            std::cout << "  normalized score:    " << mSiteScores.get_norm_score(rnaSitePosMap[i][j]);
            std::cout << std::endl << std::endl;

            ++count;
        }

    }

    return 0;

}

} // namespace rh2mfe

