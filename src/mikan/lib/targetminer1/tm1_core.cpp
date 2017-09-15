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
int TM1Core::write_site_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &mSeedTypes = mSeedSites.get_seed_types();

    seqan::CharString seedType;
    int seedStart, seedEnd;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mSeedSites.mEffectiveSites[i]) {
            continue;
        }

        seedType = mSeedTypes[i];
        seedStart = sitePos[i];
        seedEnd = seedStart + 6;

        mOFile1 << toCString(pMiRNAId) << "\t";
        mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile1 << seedStart << "\t";
        mOFile1 << seedEnd << "\t";
        mOFile1 << mSeedTypes[i] << "\t";
        mOFile1 << 0;
        mOFile1 << std::endl;

    }

    return 0;

}

int TM1Core::write_rna_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &scores = mRNAScores.get_scores();
    const seqan::String<int> &predictions = mRNAScores.get_labels();
    const seqan::String<unsigned> &siteNum = mRNAScores.get_site_num();
    mikan::TMRNAPosSet const &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();
    typedef std::multimap<double, unsigned>::iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;
    TItMap itPos;
    std::multimap<double, unsigned> sortedMRNAByScore;

    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); ++i) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }
        sortedMRNAByScore.insert(TPosPair(-1.0f * scores[i], i));
    }

    for (itPos = sortedMRNAByScore.begin(); itPos != sortedMRNAByScore.end(); ++itPos) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) mMRNAIds[uniqRNAPosSet[(*itPos).second]]) << "\t";
        mOFile2 << scores[(*itPos).second] << "\t";
        mOFile2 << siteNum[(*itPos).second] << "\t";
        mOFile2 << predictions[(*itPos).second];
        mOFile2 << std::endl;
    }

    return 0;
}

int TM1Core::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const mikan::TCharSet &mSeedTypes = mSeedSites.get_seed_types();

    seqan::CharString seedType;
    int seedStart, seedEnd;
    int count = 0;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mSeedSites.mEffectiveSites[i]) {
            continue;
        }

        seedType = mSeedTypes[i];
        seedStart = sitePos[i];
        seedEnd = seedStart + 6;

        std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
        mSiteScores.write_alignment(i);
        std::cout << "  miRNA:                " << toCString(pMiRNAId) << std::endl;
        std::cout << "  mRNA:                 " << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << std::endl;
        std::cout << "  seed type:            " << toCString(seedType) << std::endl;
        std::cout << std::endl;
        std::cout << "  position(seed start): " << seedStart << std::endl;
        std::cout << "  position(seed end):   " << seedEnd << std::endl;
        std::cout << std::endl << std::endl;

        ++count;

    }

    return 0;
}

} // namespace tm1p
