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
#include "tm1_option.hpp"        // TM1CSOptions
#include "tm1_core.hpp"          // TM1Core

namespace tm1p {

int TM1CoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line
    tm1p::TM1CSOptions options;
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
    tm1p::TM1Core tm1Core(options, mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    tm1Core.open_output_file();
    retVal = tm1Core.calculate_all_scores();
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

//
// TM1Core methods
//
void TM1Core::init_from_args(mikan::MKOptions const &opts) {
    mOutputAlign = opts.mOutputAlign;
    mOFileSite = opts.mOFileSite;
    mOFileScore = opts.mOFileTotal;

    resize(mSeedTypeDef, 1);
    mSeedTypeDef[0] = "";
}

int TM1Core::open_output_file() {
    // Open output file 1
    mOFile1.open(toCString(mOFileSite), std::ofstream::out);
    if (!mOFile1.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileSite) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    // Open output file 2
    mOFile2.open(toCString(mOFileScore), std::ofstream::out);
    if (!mOFile2.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileScore) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    return 0;
}

int TM1Core::calculate_all_scores() {
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

int TM1Core::calculate_mirna_scores(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    // Generate seed sequences
    TM1SeedSeqs seedSeqs;
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

    // Filter overlapped sites
    mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
    if (mExecFilterOverlap) {
        retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Check overlapped sites failed." << std::endl;
            return 1;
        }
    }

    // Get raw features
    if (mExecCalcSiteScore) {
        retVal = mSiteScores.calc_scores(miRNASeq, mMRNASeqs, mSeedSites, mRNAWithSites);
        if (retVal != 0) {
            std::cerr << "ERROR: Feature calculation failed." << std::endl;
            return 1;
        }
    }


    // Summarize classified results
    if (mExecSumScores) {
        retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Summarizing classified results failed." << std::endl;
            return 1;
        }
    }

    // Write site positions
    if (mOutputSitePos) {
        retVal = write_site_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write site positions." << std::endl;
            return 1;
        }
    }

    // Write total scores
    if (mOutputScore) {
        retVal = write_rna_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write scores." << std::endl;
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
    mSiteScores.clear_scores();
    mRNAWithSites.clear_maps();
    mRNAScores.clear_scores();

    return 0;
}

int TM1Core::write_site_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &mSeedTypes = mSeedSites.get_seed_types();

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
    const seqan::StringSet<seqan::CharString> &mSeedTypes = mSeedSites.get_seed_types();

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
