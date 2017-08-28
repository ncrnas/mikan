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
void MKECore::set_effective_tools() {
    resize(mEffectiveTools, mkens::TOOL_NUM, false);

    unsigned n = 0;

    mEffectiveTools[0] = true;
    mIdxMap[0] = n;
    ++n;

    mEffectiveTools[1] = true;
    mIdxMap[1] = n;
    ++n;

    mEffectiveTools[2] = true;
    mIdxMap[2] = n;
    ++n;

    mEffectiveTools[3] = true;
    mIdxMap[3] = n;
    ++n;

    mEffectiveTools[4] = true;
    mIdxMap[4] = n;
    ++n;

    mEffectiveTools[5] = true;
    mIdxMap[5] = n;
    ++n;

    mEffectiveToolN = n;
    mSeedSites.init_site_list(mEffectiveToolN);

}

void MKECore::init_score_lists() {
    for (unsigned i = 0; i < mkens::TOOL_NUM; i++) {
        if (!mEffectiveTools[i]) {
            continue;
        }
        mikan::MKCoreBase &core = get_tool_core(i);
        mikan::MKSiteScores &siteScores = core.get_site_scores();
        seqan::CharString prefix = mMKEOpts.mToolPrefix[i];
        mSiteScores.add_score_types(mMKEOpts, siteScores, prefix);
    }

}

void MKECore::clear_all() {
    for (unsigned i = 0; i < mkens::TOOL_NUM; i++) {
        mikan::MKCoreBase &core = get_tool_core(i);
        core.clear_all();
    }
    mSeedSites.clear_pos();
    mSiteScores.clear_scores();
    mRNAScores.clear_scores();
}


int MKECore::find_seed_sites(unsigned pIdx) {
    if (mFindSeedSites) {
        combine_site_pos(pIdx);
        combine_seed_types();
    }

    return 0;
}

int MKECore::combine_site_pos(unsigned pIdx) {
    for (unsigned i = 0; i < mkens::TOOL_NUM; i++) {
        if (!mEffectiveTools[i]) {
            continue;
        }

        mikan::MKCoreBase &core = get_tool_core(i);
        int retVal = core.find_seed_sites(pIdx);
        if (retVal != 0) {
            return 1;
        }

        seqan::CharString prefix = mMKEOpts.mToolPrefix[i];
        mikan::MKSeedSites &seedSites = core.get_seed_sites();
        mSeedSites.add_to_set(seedSites, mIdxMap[i], prefix);
    }

    mSeedSites.create_pos_map();

    return 0;
}

int MKECore::combine_seed_types() {
    for (unsigned i = 0; i < mkens::TOOL_NUM; i++) {
        if (!mEffectiveTools[i]) {
            continue;
        }

        mikan::MKCoreBase &core = get_tool_core(i);
        mikan::MKSeedSites &seedSites = core.get_seed_sites();
        seqan::CharString prefix = mMKEOpts.mToolPrefix[i];
        mSeedSites.add_seed_types(seedSites, mIdxMap[i], prefix);
    }

    mSeedSites.combine_seed_types();

    return 0;
}


int MKECore::calc_site_scores(unsigned pIdx) {
    if (mCalcSiteScore) {
        mSiteScores.init_score_list(mSeedSites);

        for (unsigned i = 0; i < mkens::TOOL_NUM; i++) {
            if (!mEffectiveTools[i]) {
                continue;
            }

            mikan::MKCoreBase &core = get_tool_core(i);
            int retVal = core.calc_site_scores(pIdx);
            if (retVal != 0) {
                return 1;
            }

            mikan::MKSeedSites &seedSites = core.get_seed_sites();
            mikan::MKSiteScores &seedScores = core.get_site_scores();
            seqan::CharString prefix = mMKEOpts.mToolPrefix[i];

            mSiteScores.add_scores(mMKEOpts, seedSites, mSeedSites, seedScores, prefix);

        }

        mSiteScores.combine_scores(mMKEOpts);

    }

    return 0;
}

int MKECore::calc_rna_scores(unsigned pIdx) {
    int retVal;

    if (mCalcRNAScore) {
        for (unsigned i = 0; i < mkens::TOOL_NUM; i++) {
            if (!mEffectiveTools[i]) {
                continue;
            }

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

int MKECore::write_site_score(seqan::CharString const &pMiRNAId) {

    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    seqan::CharString seedType;
    unsigned posA1;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mSeedSites.mEffectiveSites[i] || !mSiteScores.mEffectiveSites[i]) {
            continue;
        }

        posA1 = sitePos[i];

        mOFile1 << toCString(pMiRNAId) << "\t";
        mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile1 << posA1 << "\t";
        mOFile1 << seedTypes[i] << "\t";
        mOFile1 << mSiteScores.get_score(i) << "\t";
        mOFile1 << toCString(mSiteScores.get_tool_score(i));
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

