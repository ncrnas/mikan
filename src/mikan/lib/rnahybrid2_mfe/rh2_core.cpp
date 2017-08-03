#include <math.h>                // roundf
#include <iostream>
#include <string>                // string
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                 // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"          // MKInput
#include "rh2_option.hpp"        // RH2Options
#include "rh2_seed_site.hpp"     // RH2Sequences, RH2SeedSites
#include "rh2_site_score.hpp"    // RH2SiteScores, RH2TotalScores
#include "rh2_site_filter.hpp"   // RH2Overlap, RH2TopNSites, RH2SortedSitePos
#include "rh2_core.hpp"          // RH2Core

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
    mikan::MKInput coreInput;
    coreInput.set_options(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0) {
        return 1;
    }

    // Create index
    mikan::TRNASet const &mMRNASeqs = coreInput.get_mrna_seqs();
    mikan::TIndexQGram index(mMRNASeqs);
    mikan::TFinder finder(index);

    // Calculate scores for all miRNAs
    mikan::TCharSet const &mMiRNAIds = coreInput.get_mirna_ids();
    mikan::TRNASet const &mMiRNASeqs = coreInput.get_mirna_seqs();
    mikan::TCharSet const &mMRNAIds = coreInput.get_mrna_ids();
    rh2mfe::RH2Core rh2Core(options, mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
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
void RH2Core::init_from_args(mikan::MKOptions const &opts) {
    mOutputAlign = opts.mOutputAlign;
    mOFileMFE = opts.mOFileSite;
    mOFileTotal = opts.mOFileTotal;

    resize(mSeedTypeDef, 1);
    mSeedTypeDef[0] = opts.mSeedDef;
    mMaxHits = opts.mMaxHits;

    mOverlapDef = opts.mOverlapDef;

    mSiteScores.init_from_args();
    mSiteFilter.init_from_args();
    mTopNSites.init_from_args();

}

int RH2Core::open_output_file() {
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

int RH2Core::calculate_all_scores() {
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

int RH2Core::calculate_mirna_scores(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    // Generate seed sequences
    RH2SeedSeqs seedSeqs;
    seedSeqs.set_mirna_seq(miRNASeq);
    seedSeqs.set_flags(mSeedTypeDef);
    retVal = seedSeqs.create_seed_seqs();
    if (retVal != 0) {
        std::cerr << "ERROR: Generate seed sequences failed." << std::endl;
        return 1;
    }

    // Search seed sites
    if (mExecSearchSeedSites) {
        retVal = mSeedSites.find_seed_sites(seedSeqs, mSeedTypeDef);
        if (retVal != 0) {
            std::cerr << "ERROR: Seed site search failed." << std::endl;
            return 1;
        }
    }

    // Calculate MFE values
    if (mExecCalSiteScore) {
        retVal = mSiteScores.calc_scores(miRNASeq, mMRNASeqs, mSeedSites, mRNAWithSites);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate MFE values failed." << std::endl;
            return 1;
        }
    }

    // Filter overlapped sites
    mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
    if (mExecFilterOverlap) {
        retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Check overlapped sites failed." << std::endl;
            return 1;
        }
    }

    // Filter sites by the numbers of sites
    if (mExecFilterSiteNum) {
        retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Filter top scored sites failed." << std::endl;
            return 1;
        }
    }

    // Summarize MFE values
    if (mExecSumScores) {
        retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
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
    mRNAWithSites.clear_maps();
    mSiteScores.clear_scores();
    mRNAScores.clear_scores();

    return 0;
}

int RH2Core::write_mfe_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    for (unsigned i = 0; i < length(mRNAWithSites.mEffectiveRNAs); i++) {
        if (!mRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        int seedStart;
        float score;
        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!mSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }

            seedStart = sitePos[rnaSitePosMap[i][j]];
            score = mSiteScores.get_score(rnaSitePosMap[i][j]);
            score = roundf(score * 10.0f) / 10.0f;

            mOFile1 << toCString(pMiRNAId) << "\t";
            mOFile1 << toCString((seqan::CharString) mMRNAIds[uniqRNAPosSet[i]]) << "\t";
            mOFile1 << seedStart + 1 << "\t";
            mOFile1 << seedStart + 7 << "\t";
            //        mOFile1 << mSiteScores.get_wide_site_start(posIdx) + 1  << "\t";
            mOFile1 << toCString((seqan::CharString) seedTypes[rnaSitePosMap[i][j]]) << "\t";
            mOFile1 << score << "\t";
            mOFile1 << mSiteScores.get_norm_score(rnaSitePosMap[i][j]);
            mOFile1 << std::endl;
        }

    }

    return 0;

}

int RH2Core::write_total_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &totalScores = mRNAScores.get_scores();
    const seqan::String<float> &totalNormScores = mRNAScores.get_norm_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();

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

int RH2Core::write_alignment(seqan::CharString const &pMiRNAId) {
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
            score = roundf(score * 10.0f) / 10.0f;

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.write_alignment(rnaSitePosMap[i][j]);
            std::cout << "  miRNA:               " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:                " << toCString((seqan::CharString) mMRNAIds[uniqRNAPosSet[i]])
                      << std::endl;
            std::cout << "  seed type:           " << toCString((seqan::CharString) seedTypes[rnaSitePosMap[i][j]])
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

