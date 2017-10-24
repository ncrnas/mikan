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
void TSSVMCore::write_site_score_tab(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    if (mPrintSiteHeader) {
        mOFile1 << "# miRNA name, ";
        mOFile1 << "mRNA name, ";
        mOFile1 << "start (1-base), ";
        mOFile1 << "end (1-base), ";
        mOFile1 << "seed type, ";
        mOFile1 << "score 1 (discriminant value), ";
        mOFile1 << "score 2 (not used)";
        mOFile1 << std::endl;
        mPrintSiteHeader = false;
    }

    mikan::TCharStr seedType = seedTypes[pSitePosIdx];
    int seedStart = sitePos[pSitePosIdx];
    if (seedType == "7mer-A1") {
        seedStart += 1;
    }

    float score = mSiteScores.get_score(pSitePosIdx);
    score = roundf(score * 10000.0f) / 10000.0f;
    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart + 1 << "\t";
    mOFile1 << seedStart + 7 << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << score << "\t";
    mOFile1 << 0;
    mOFile1 << std::endl;

}

void TSSVMCore::write_site_score_gff(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();

    mikan::TCharStr seedType = seedTypes[pSitePosIdx];
    int seedStart = sitePos[pSitePosIdx];
    if (seedType == "7mer-A1") {
        seedStart += 1;
    }

    float score = mSiteScores.get_score(pSitePosIdx);
    score = roundf(score * 10000.0f) / 10000.0f;
    mOFile1 << toCString(pMiRNAId) << "\t";
    mOFile1 << toCString((mikan::TCharStr) mMRNAIds[pRNAPosIdx]) << "\t";
    mOFile1 << seedStart + 1 << "\t";
    mOFile1 << seedStart + 7 << "\t";
    mOFile1 << toCString((mikan::TCharStr) seedTypes[pSitePosIdx]) << "\t";
    mOFile1 << score << "\t";
    mOFile1 << std::endl;

}

void TSSVMCore::write_rna_score_tab(mikan::TCharStr const &pMiRNAId) {
    typedef std::multimap<float, unsigned>::reverse_iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;

    const seqan::String<float> &scores = mRNAScores.get_scores();
    const seqan::String<int> &siteCount = mRNAScores.get_site_num();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    if (mPrintRNAheader) {
        mOFile2 << "# miRNA name, ";
        mOFile2 << "mRNA name, ";
        mOFile2 << "number of sites, ";
        mOFile2 << "score 1 (discriminant value), ";
        mOFile2 << "score 2 (not used)";
        mOFile2 << std::endl;
        mPrintRNAheader = false;
    }

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
        mOFile2 << toCString((mikan::TCharStr) mMRNAIds[uniqRNAPosSet[(*itPos).second]]) << "\t";
        mOFile2 << siteCount[(*itPos).second] << "\t";
        mOFile2 << score << "\t";
        mOFile2 << 0;
        mOFile2 << std::endl;
    }
}

int TSSVMCore::write_alignment(mikan::TCharStr const &pMiRNAId) {

    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<float> &scores = mSiteScores.get_scores();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();

    unsigned padw = 19;
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
            if (seedType == "7mer-A1") {
                seedStart += 1;
            }

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.write_alignment(rnaSitePosMap[i][j]);
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "miRNA: " << toCString(pMiRNAId) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "mRNA: ";
            std::cout << toCString((mikan::TCharStr) mMRNAIds[mRNAPos[rnaSitePosMap[i][j]]]) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "seed type: " << toCString(seedType) << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "start (1-base): " << seedStart + 1 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "end (1-base): " << seedStart + 7 << std::endl;
            std::cout << std::right << std::setw(padw) << std::setfill(' ');
            std::cout << "site level score: " << scores[rnaSitePosMap[i][j]] << std::endl;
            std::cout << std::endl;

            ++count;
        }
    }

    return 0;

}

} // namespace tssvm
