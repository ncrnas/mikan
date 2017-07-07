#include <math.h>                // roundf
#include <iostream>
#include <string>                // string
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                 // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_inst_template.hpp"  // TRNATYPE
#include "rh2_option.hpp"        // RH2Options
#include "rh2_seed_site.hpp"     // RH2Sequences, RH2SeedSites
#include "rh2_score.hpp"         // RH2MFEScores, RH2TotalScores
#include "rh2_site_cluster.hpp"  // RH2Overlap, RH2TopNScore, RH2SortedSitePos
#include "rh2_core.hpp"          // RH2Core
#include "mk_input.hpp"          // MKInput

using namespace mikan;

namespace rh2mfe {

int RH2CoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line.
    rh2mfe::RH2Options options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK) {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    mikan::MKInput<TRNATYPE> coreInput;
    coreInput.set_options(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0) {
        return 1;
    }

    // Create index
    rh2mfe::RH2Core<TRNATYPE>::TRNASet const &mMRNASeqs = coreInput.get_mrna_seqs();
    rh2mfe::RH2Core<TRNATYPE>::TIndexQGram index(mMRNASeqs);
    rh2mfe::RH2Core<TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    rh2mfe::RH2Core<TRNATYPE>::TCharSet const &mMiRNAIds = coreInput.get_mirna_ids();
    rh2mfe::RH2Core<TRNATYPE>::TRNASet const &mMiRNASeqs = coreInput.get_mirna_seqs();
    rh2mfe::RH2Core<TRNATYPE>::TCharSet const &mMRNAIds = coreInput.get_mrna_ids();
    int mRNAMaxLen = options.mTargetLen;
    int miRNAMaxLen = options.mQueryLen;
    std::string seedDef(toCString(options.mSeedDef));
    rh2mfe::RH2Core<TRNATYPE> rh2Core(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder,
                                      mRNAMaxLen, miRNAMaxLen, seedDef);
    rh2Core.init_from_args(options);
    rh2Core.open_output_file();
    retVal = rh2Core.calculate_all_scores();
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

//
// RH2Core methods
//
template<class TRNAString, int SEEDLEN>
void RH2Core<TRNAString, SEEDLEN>::init_from_args(RH2Options &opts) {
    mOutputAlign = opts.mOutputAlign;
    mOFileMFE = opts.mOFileSite;
    mOFileTotal = opts.mOFileTotal;

    mSeedDef = opts.mSeedDef;
    mOverlapDef = opts.mOverlapDef;
    mMaxHits = opts.mMaxHits;

}

template<class TRNAString, int SEEDLEN>
int RH2Core<TRNAString, SEEDLEN>::open_output_file() {
    // Open output file 1
    mOFile1.open(toCString(mOFileMFE), std::ofstream::out);
    if (!mOFile1.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileMFE) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    // Open output file 2
    mOFile2.open(toCString(mOFileTotal), std::ofstream::out);
    if (!mOFile2.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileTotal) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int RH2Core<TRNAString, SEEDLEN>::calculate_all_scores() {
    int retVal;

    for (unsigned i = 0; i < length(mMiRNASeqs); ++i) {

#if SEQAN_ENABLE_DEBUG
        clock_t startTime = clock();
#endif

        retVal = calculate_mirna_scores(i);
        if (retVal != 0) {
            std::cerr << "ERROR: Score calculation failed for " << toCString((seqan::CharString) mMiRNAIds[i]);
            std::cerr << "." << std::endl;
            return 1;
        }

#if SEQAN_ENABLE_DEBUG
        std::cout << toCString((seqan::CharString) mMiRNAIds[i]) << ": ";
        std::cout << double(clock() - startTime) / (double) CLOCKS_PER_SEC << " seconds." << std::endl;
#endif

    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int RH2Core<TRNAString, SEEDLEN>::calculate_mirna_scores(unsigned pIdx) {
    int retVal;
    TRNAString miRNASeq = mMiRNASeqs[pIdx];

    // Search seed sites
    if (mExecSearchSeedSites) {
        retVal = mSeedSites.find_seed_sites(miRNASeq, mSeedDef, mOverlapDef);
        if (retVal != 0) {
            std::cerr << "ERROR: Seed site search failed." << std::endl;
            return 1;
        }
    }

    // Calculate MFE values
    if (mExecCalMFEScore) {
        retVal = mMfeScores.calc_scores(mSeedSites, miRNASeq, mMRNASeqs, mOverlapDef);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate MFE values failed." << std::endl;
            return 1;
        }
    }

    // Filter overlapped sites
    if (mExecFilterOverlap) {
        retVal = mOverlappedSites.filter_overlapped_sites(mSeedSites, mMfeScores, mOverlapDef);
        if (retVal != 0) {
            std::cerr << "ERROR: Check overlapped sites failed." << std::endl;
            return 1;
        }
    }

    // Filter sites by the numbers of sites
    if (mExecFilterSiteNum) {
        retVal = mTopScoredSites.filter_sites(mSeedSites, mMfeScores, mMaxHits);
        if (retVal != 0) {
            std::cerr << "ERROR: Filter top scored sites failed." << std::endl;
            return 1;
        }
    }

    // Sort target sites
    if (mExecSortSites) {
        retVal = mSortedSites.generate_sorted_mrna_pos(mSeedSites, mMfeScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Sort target sites failed." << std::endl;
            return 1;
        }
    }

    // Summarize MFE values
    if (mExecSumScores) {
        const seqan::String<unsigned> &sortedPos = mSortedSites.get_sorted_mrna_pos();
        retVal = mTotalScores.calc_scores(mSeedSites, mMfeScores, sortedPos);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate total MFE values failed." << std::endl;
            return 1;
        }
    }

    // Write MFE values
    if (mOutputMFEScore) {
        retVal = write_mfe_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write MFE scores." << std::endl;
            return 1;
        }
    }

    // Write total MEF values
    if (mOutputTotalScore) {
        retVal = write_total_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write total MFE scores." << std::endl;
            return 1;
        }
    }

    // Write alignments
    if (mOutputAlign) {
        retVal = write_alignment(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write alignments." << std::endl;
            return 1;
        }
    }

    mSeedSites.clear_pos();
    mMfeScores.clear_scores();
    mOverlappedSites.clear_cluster();
    mTopScoredSites.clear_cluster();
    mSortedSites.clear_site_pos();
    mTotalScores.clear_scores();

    return 0;
}

template<class TRNAString, int SEEDLEN>
int RH2Core<TRNAString, SEEDLEN>::write_mfe_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sortedPos = mSortedSites.get_sorted_mrna_pos();
    seqan::CharString seedType;
    int seedStart = 0;
    int posIdx;
    float score;

    for (unsigned i = 0; i < length(sortedPos); ++i) {
        posIdx = sortedPos[i];
        if (!mMfeScores.mEffectiveSites[posIdx]) {
            continue;
        }

        seedStart = sitePos[posIdx];
        score = mMfeScores.get_score(posIdx);
        score = roundf(score * 10.0f) / 10.0f;

        mOFile1 << toCString(pMiRNAId) << "\t";
        mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[posIdx]]) << "\t";
        mOFile1 << seedStart + 1 << "\t";
        mOFile1 << seedStart + 7 << "\t";
        //        mOFile1 << mMfeScores.get_hit_start(posIdx) + 1  << "\t";
        mOFile1 << toCString((seqan::CharString) seedTypes[posIdx]) << "\t";
        mOFile1 << score << "\t";
        mOFile1 << mMfeScores.get_norm_score(posIdx);
        mOFile1 << std::endl;
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int RH2Core<TRNAString, SEEDLEN>::write_total_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &totalScores = mTotalScores.get_scores();
    const seqan::String<float> &totalNormScores = mTotalScores.get_norm_scores();
    const seqan::String<int> &mRNAPos = mTotalScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mTotalScores.get_site_num();

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile2 << totalScores[i] << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << totalNormScores[i] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int RH2Core<TRNAString, SEEDLEN>::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<unsigned> &sortedPos = mSortedSites.get_sorted_mrna_pos();
    seqan::CharString seedType;
    int seedStart = 0;
    int posIdx;
    float score;
    int count = 0;

    for (unsigned i = 0; i < length(sortedPos); ++i) {
        posIdx = sortedPos[i];

        if (!mMfeScores.mEffectiveSites[posIdx]) {
            continue;
        }

        seedStart = sitePos[posIdx];
        score = mMfeScores.get_score(posIdx);
        score = roundf(score * 10.0f) / 10.0f;

        std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
        mMfeScores.write_alignment(posIdx);
        std::cout << "  miRNA:               " << toCString(pMiRNAId) << std::endl;
        std::cout << "  mRNA:                " << toCString((seqan::CharString) mMRNAIds[mRNAPos[posIdx]])
                  << std::endl;
        std::cout << "  seed type:           " << toCString((seqan::CharString) seedTypes[posIdx]) << std::endl;
        std::cout << "  position(target 5'): " << mMfeScores.get_hit_start(posIdx) + 1;
        std::cout << std::endl;
        std::cout << "  position(seed):      " << seedStart + 1 << std::endl;
        std::cout << "  mfe:                 " << score << " kcal/mol" << std::endl;
        std::cout << "  normalized score:    " << mMfeScores.get_norm_score(posIdx);
        std::cout << std::endl << std::endl;

        ++count;

    }

    return 0;
}

// Explicit template instantiation

template
class RH2Core<TRNATYPE>;

} // namespace rh2mfe

