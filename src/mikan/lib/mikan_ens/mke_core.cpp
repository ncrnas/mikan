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
#include "mke_seed_site.hpp"     // MKESeedSites
#include "mke_site_score.hpp"    // MKESiteScores
#include "mke_core.hpp"          // MKECore

namespace mkens {

//
// MKECore methods
//
int MKECore::find_seed_sites(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    if (mFindSeedSites) {
        for (unsigned i = 0; i < TOOL_NUM; i++) {
            mikan::MKCoreBase &core = get_tool_core(i);
            retVal = core.find_seed_sites(pIdx);
            if (retVal != 0) {
                return 1;
            }
        }
    }

    return 0;
}

int MKECore::calc_site_scores(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    if (mCalcSiteScore) {
        for (unsigned i = 0; i < TOOL_NUM; i++) {
            mikan::MKCoreBase &core = get_tool_core(i);
            retVal = core.calc_site_scores(pIdx);
            if (retVal != 0) {
                return 1;
            }
        }
    }

    return 0;
}

int MKECore::calc_rna_scores(unsigned pIdx) {
    int retVal;

    if (mCalcRNAScore) {
        for (unsigned i = 0; i < TOOL_NUM; i++) {
            mikan::MKCoreBase &core = get_tool_core(i);
            retVal = core.calc_rna_scores(pIdx);
            if (retVal != 0) {
                return 1;
            }
        }
    }

    return 0;
}

int MKECore::output_results(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    // Write site scores
    if (mOutputSite) {
        retVal = write_site_score(mMiRNAIds[pIdx]);
        if (retVal != 0) { ;
            return 1;
        }
    }

    // Write total scores
    if (mOutputRNA) {
        retVal = write_rna_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;
}

void MKECore::clear_all() {
    for (unsigned i = 0; i < TOOL_NUM; i++) {
        mikan::MKCoreBase &core = get_tool_core(i);
        core.clear_all();
    }
    mSeedSites.clear_pos();
    mSiteScores.clear_scores();
    mRNAScores.clear_scores();
}

int MKECore::write_site_score(seqan::CharString const &pMiRNAId) {

    const mikan::TCharSet &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    seqan::CharString seedType;
    int seedStart, seedEnd;


    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mSiteScores.mEffectiveSites[i]) {
            continue;
        }

        seedStart = sitePos[i];
        seedEnd = seedStart + 6;

        mOFile1 << toCString(pMiRNAId) << "\t";
        mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile1 << seedStart << "\t";
        mOFile1 << seedEnd << "\t";
        mOFile1 << seedTypes[i] << "\t";
        mOFile1 << mSiteScores.get_score(i);
        mOFile1 << std::endl;
    }

    return 0;
}

int MKECore::write_rna_score(seqan::CharString const &pMiRNAId) {

    const seqan::String<float> &totalScores = mRNAScores.get_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();
    typedef std::multimap<double, unsigned>::iterator TItMap;
    typedef std::pair<double, unsigned> TPosPair;
    TItMap itPos;
    std::multimap<double, unsigned> sortedMRNAByScore;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        sortedMRNAByScore.insert(TPosPair((float) totalScores[i], i));
    }

    for (itPos = sortedMRNAByScore.begin(); itPos != sortedMRNAByScore.end(); ++itPos) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) mMRNAIds[mRNAPos[(*itPos).second]]) << "\t";
        mOFile2 << totalScores[(*itPos).second] << "\t";
        mOFile2 << siteNum[(*itPos).second];
        mOFile2 << std::endl;
    }

    return 0;
}


} // namespace mkens

