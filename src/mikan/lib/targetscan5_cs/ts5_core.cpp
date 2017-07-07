#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                 // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <map>                   // multimap
#include <utility>               // pair
#include <seqan/arg_parse.h>
#include <ts5_inst_template.hpp> // TRNATYPE
#include <ts5_option.hpp>        // TS5CSOptions
#include <ts5_seed_site.hpp>     // TS5SeedSites
#include <ts5_feature.hpp>       // TS5RawFeatures
#include <ts5_score.hpp>         // TS5ContextScores, TS5TotalScores
#include <ts5_core.hpp>          // TS5Core
#include <mk_input.hpp>          // MKInput

namespace ts5cs {

int TS5CoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line
    ts5cs::TS5CSOptions options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK) {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    mikan::MKInput<ts5cs::TRNATYPE> coreInput;
    coreInput.set_file_names(options.mMiRNAFasta, options.mMRNAFasta);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0) {
        return 1;
    }

    // Create index
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TRNASet const &mMRNASeqs = coreInput.get_mrna_seqs();
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TCharSet const &mMiRNAIds = coreInput.get_mirna_ids();
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TRNASet const &mMiRNASeqs = coreInput.get_mirna_seqs();
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TCharSet const &mMRNAIds = coreInput.get_mrna_ids();
    ts5cs::TS5Core<ts5cs::TRNATYPE> ts5Core(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    ts5Core.init_from_args(options);
    ts5Core.open_output_file();
    retVal = ts5Core.calculate_all_scores();
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

//
// TS5Core methods
//
template<class TRNAString, int SEEDLEN>
void TS5Core<TRNAString, SEEDLEN>::init_from_args(TS5CSOptions &opts) {
    mOutputAlign = opts.mOutputAlign;
    mOFileContext = opts.mOFileSite;
    mOFileTotal = opts.mOFileTotal;
}

template<class TRNAString, int SEEDLEN>
int TS5Core<TRNAString, SEEDLEN>::open_output_file() {
    // Open output file 1
    mOFile1.open(toCString(mOFileContext), std::ofstream::out);
    if (!mOFile1.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileContext) << std::endl;
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
int TS5Core<TRNAString, SEEDLEN>::calculate_all_scores() {
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
int TS5Core<TRNAString, SEEDLEN>::calculate_mirna_scores(unsigned pIdx) {
    int retVal;

    // Search seed sites
    if (mExecSearchSeedSites) {
        retVal = mSeedSites.find_seed_sites(mMiRNASeqs[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Seed site search failed." << std::endl;
            return 1;
        }
    }

    // Get raw features
    if (mExecGetRawFeat) {
        retVal = mRawFeatures.add_features(mMiRNASeqs[pIdx], mMRNASeqs, mSeedSites);
        if (retVal != 0) {
            std::cerr << "ERROR: Feature calculation failed." << std::endl;
            return 1;
        }
    }

    // Calculate context scores
    if (mExecCalcContexScore) {
        retVal = mCsScores.calc_scores(mRawFeatures);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate regression values failed." << std::endl;
            return 1;
        }
    }

    // Summarize context scores
    if (mExecSumScores) {
        retVal = mTotalScore.calc_scores(mSeedSites, mCsScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate total regression values failed." << std::endl;
            return 1;
        }
    }

    // Write context scores
    if (mOutputContexScore) {
        retVal = write_context_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write context scores." << std::endl;
            return 1;
        }
    }

    // Write total scores
    if (mOutputTotalScore) {
        retVal = write_total_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write total scores." << std::endl;
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
    mRawFeatures.clear_features();
    mCsScores.clear_scores();
    mTotalScore.clear_scores();

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TS5Core<TRNAString, SEEDLEN>::write_context_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    seqan::CharString seedType;
    int seedStart, seedEnd;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mCsScores.mEffectiveSites[i]) {
            continue;
        }

        seedType = mRawFeatures.get_seed_type(i);
        seedStart = sitePos[i];
        if (seedType == "7mer-A1") {
            seedStart += 1;
        }

        seedEnd = seedStart + 6;
        if (seedType == "8mer") {
            seedEnd += 1;
        }

        mOFile1 << toCString(pMiRNAId) << "\t";
        mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << "\t";
        mOFile1 << seedStart << "\t";
        mOFile1 << seedEnd << "\t";
        mOFile1 << toCString(mRawFeatures.get_seed_type(i)) << "\t";
        mOFile1 << mCsScores.get_score(i);
        mOFile1 << std::endl;
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TS5Core<TRNAString, SEEDLEN>::write_total_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &totalScores = mTotalScore.get_scores();
    const seqan::String<int> &mRNAPos = mTotalScore.get_mrna_pos();
    const seqan::String<int> &siteNum = mTotalScore.get_site_num();
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

template<class TRNAString, int SEEDLEN>
int TS5Core<TRNAString, SEEDLEN>::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const TS5Alignment<TRNAString> &alignment = mRawFeatures.get_alignment();
    seqan::CharString seedType;
    int seedStart, seedEnd;
    int count = 0;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        if (!mCsScores.mEffectiveSites[i]) {
            continue;
        }

        seedType = mRawFeatures.get_seed_type(i);
        seedStart = sitePos[i];
        if (seedType == "7mer-A1") {
            seedStart += 1;
        }

        seedEnd = seedStart + 6;
        if (seedType == "8mer") {
            seedEnd += 1;
        }

        std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
        alignment.write_alignment(i);
        std::cout << "  miRNA:                " << toCString(pMiRNAId) << std::endl;
        std::cout << "  mRNA:                 " << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << std::endl;
        std::cout << "  seed type:            " << toCString(seedType) << std::endl;
        std::cout << "  position(seed start): " << seedStart << std::endl;
        std::cout << "  position(seed end):   " << seedEnd << std::endl;
        std::cout << "  context score:        " << mCsScores.get_score(i);
        std::cout << std::endl << std::endl;

        ++count;

    }

    return 0;
}

// Explicit template instantiation

template
class TS5Core<TRNATYPE>;

} // namespace ts5cs
