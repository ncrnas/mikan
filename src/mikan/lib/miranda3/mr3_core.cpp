#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"           // MKInput
#include "mr3_option.hpp"         // MR3Options
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "mr3_site_score.hpp"     // MR3SiteScores
#include "mr3_site_filter.hpp"    // MR3SiteFilter
#include "mr3_core.hpp"           // MR3Core

namespace mr3as {

int MR3CoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line.
    mr3as::MR3Options options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK) {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    mikan::MKInput coreInput;
    coreInput.set_options(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0) {
        return retVal;
    }

    // Create index
    mikan::TRNASet const &mMRNASeqs = coreInput.get_mrna_seqs();
    mikan::TIndexQGram index(mMRNASeqs);
    mikan::TFinder finder(index);

    // Calculate scores for all miRNAs
    mikan::TCharSet const &mMiRNAIds = coreInput.get_mirna_ids();
    mikan::TRNASet const &mMiRNASeqs = coreInput.get_mirna_seqs();
    mikan::TCharSet const &mMRNAIds = coreInput.get_mrna_ids();

    mr3as::MR3Core mr3Core(options, mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    mr3Core.open_output_file();
    retVal = mr3Core.calculate_all_scores();

    return retVal;
}

//
// MR3Core methods
//
void MR3Core::init_from_args(mikan::MKOptions const &opts) {
    mOutputAlign = opts.mOutputAlign;
    mOFileSite = opts.mOFileSite;
    mOFileRNA = opts.mOFileTotal;

    resize(mSeedTypeDef, 6);
    mSeedTypeDef[0] = 'Y';
    mSeedTypeDef[1] = 'Y';
    mSeedTypeDef[2] = 'Y';
    if (opts.mMinSeedLen == 7) {
        mSeedTypeDef[0] = 'N';
    } else if (opts.mMinSeedLen == 8) {
        mSeedTypeDef[0] = 'N';
        mSeedTypeDef[1] = 'N';
    }

    if (opts.mMaxSeedLen == 7) {
        mSeedTypeDef[2] = 'N';
    } else if (opts.mMaxSeedLen == 6) {
        mSeedTypeDef[2] = 'N';
        mSeedTypeDef[1] = 'N';
    }
    mSeedTypeDef[3] = opts.mAllowGUWobble;
    mSeedTypeDef[4] = opts.mAllowMismatch;
    mSeedTypeDef[5] = opts.mAllowBT;

}

int MR3Core::open_output_file() {
    // Open output file 1
    mOFile1.open(toCString(mOFileSite), std::ofstream::out);
    if (!mOFile1.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileSite) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    // Open output file 2
    mOFile2.open(toCString(mOFileRNA), std::ofstream::out);
    if (!mOFile2.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileRNA) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    return 0;
}

int MR3Core::calculate_all_scores() {
    int retVal;

    for (unsigned i = 0; i < length(mMiRNASeqs); ++i) {

#if SEQAN_ENABLE_DEBUG
        clock_t startTime = clock();
#endif

        retVal = calculate_mirna_scores(i);
        if (retVal != 0) {
            std::cerr << "ERROR: Score calculation failed for ";
            std::cerr << toCString((seqan::CharString) mMiRNAIds[i]) << "." << std::endl;
            return 1;
        }

#if SEQAN_ENABLE_DEBUG
        std::cout << toCString((seqan::CharString) mMiRNAIds[i]) << ": ";
        std::cout << double(clock() - startTime) / (double) CLOCKS_PER_SEC << " seconds." << std::endl;
#endif

    }

    return 0;
}

int MR3Core::calculate_mirna_scores(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    // Generate seed sequences
    mSeedSeqs.set_flags(mSeedTypeDef);
    retVal = mSeedSeqs.create_seed_seqs(miRNASeq);
    if (retVal != 0) {
        std::cerr << "ERROR: Generate seed sequences failed." << std::endl;
        return 1;
    }

    // Search seed sites
    if (mExecSearchSeedSites) {
        retVal = mSeedSites.find_seed_sites(mSeedSeqs, mSeedTypeDef);
        if (retVal != 0) {
            std::cerr << "ERROR: Seed site search failed." << std::endl;
            return 1;
        }
    }

    // Calculate alignment and energy scores
    if (mExecCalSiteScore) {
        retVal = mSiteScores.calc_scores(miRNASeq, mMRNASeqs, mSeedSites, mRNAWithSites);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate site scores failed." << std::endl;
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

    // Summarize alignment and energy scores
    if (mExecSumScores) {
        retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate total scores failed." << std::endl;
            return 1;
        }
    }

    // Write site scores
    if (mOutputSiteScore) {
        retVal = write_site_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write site scores." << std::endl;
            return 1;
        }
    }

    // Write total scores
    if (mOutputTotalScore) {
        retVal = write_rna_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write total ddG scores." << std::endl;
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

    mSeedSeqs.clear_seeds();
    mSeedSites.clear_pos();
    mSiteScores.clear_scores();
    mRNAWithSites.clear_maps();
    mRNAScores.clear_scores();

    return 0;
}

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
