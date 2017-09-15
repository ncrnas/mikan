#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                    // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"             // MKInput
#include "tssvm_seed_site.hpp"      // TSSVMSeedSites, TSSVMSiteFilter
#include "tssvm_align.hpp"          // TSAlign
#include "tssvm_site_svm.hpp"       // TSSVMSiteModel, TSSVMSiteInputVector
#include "tssvm_mrna_feature.hpp"   // TSSVMRNARawFeatures
#include "tssvm_core.hpp"           // TSSVMCore

namespace tssvm {

//
// TSSVMCore methods
//
int TSSVMCore::write_site_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<float> &scores = mSiteScores.get_scores();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();

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
            if (seedType == "7mer-A1") {
                seedStart += 1;
            }

            float score = scores[rnaSitePosMap[i][j]];
            score = roundf(score * 10000.0f) / 10000.0f;
            mOFile1 << toCString(pMiRNAId) << "\t";
            mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[rnaSitePosMap[i][j]]]) << "\t";
            mOFile1 << seedStart + 1 << "\t";
            mOFile1 << seedStart + 7 << "\t";
            mOFile1 << toCString((seqan::CharString) seedTypes[rnaSitePosMap[i][j]]) << "\t";
            mOFile1 << score << "\t";
            mOFile1 << std::endl;
        }

    }

    return 0;

}

int TSSVMCore::write_rna_score(seqan::CharString const &pMiRNAId) {
    typedef std::multimap<float, unsigned>::reverse_iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;

    const seqan::String<float> &scores = mRNAScores.get_scores();
    const seqan::String<int> &siteCount = mRNAScores.get_site_num();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    std::multimap<float, unsigned> sortedMRNAByScore;
    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); i++) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        sortedMRNAByScore.insert(TPosPair((float) scores[i], i));
    }

    for (TItMap itPos = sortedMRNAByScore.rbegin(); itPos != sortedMRNAByScore.rend(); ++itPos) {
        float score = scores[(*itPos).second];
        score = roundf(score * 10000.0f) / 10000.0f;

        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) mMRNAIds[uniqRNAPosSet[(*itPos).second]]) << "\t";
        mOFile2 << score << "\t";
        mOFile2 << siteCount[(*itPos).second] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

int TSSVMCore::write_alignment(seqan::CharString const &pMiRNAId) {

    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<float> &scores = mSiteScores.get_scores();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();

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
            if (seedType == "7mer-A1") {
                seedStart += 1;
            }

            std::cout << "### " << (i + j) + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.write_alignment(rnaSitePosMap[i][j]);
            std::cout << "  miRNA:                " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:                 ";
            std::cout << toCString((seqan::CharString) mMRNAIds[mRNAPos[rnaSitePosMap[i][j]]]) << std::endl;
            std::cout << "  seed type:            " << toCString(seedType) << std::endl;
            std::cout << "  position(seed start): " << seedStart + 1 << std::endl;
            std::cout << "  site level score:     " << scores[rnaSitePosMap[i][j]];
            std::cout << std::endl << std::endl;


        }
    }

    return 0;

}

} // namespace tssvm
