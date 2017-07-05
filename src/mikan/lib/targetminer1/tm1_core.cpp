#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                 // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <map>                   // multimap
#include <utility>               // pair
#include <seqan/arg_parse.h>
#include <tm1_inst_template.hpp> // TRNATYPE
#include <tm1_option.hpp>        // TM1CSOptions
#include <tm1_core.hpp>          // TM1CoreInput, TM1Core

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
    tm1p::TM1CoreInput<tm1p::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0) {
        return 1;
    }

    // Create index
    tm1p::TM1Core<tm1p::TRNATYPE>::TRNASet const &mMRNASeqs = coreInput.get_mrna_seqs();
    tm1p::TM1Core<tm1p::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    tm1p::TM1Core<tm1p::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    tm1p::TM1Core<tm1p::TRNATYPE>::TCharSet const &mMiRNAIds = coreInput.get_mirna_ids();
    tm1p::TM1Core<tm1p::TRNATYPE>::TRNASet const &mMiRNASeqs = coreInput.get_mirna_seqs();
    tm1p::TM1Core<tm1p::TRNATYPE>::TCharSet const &mMRNAIds = coreInput.get_mrna_ids();
    tm1p::TM1Core<tm1p::TRNATYPE> tm1Core(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    tm1Core.init_from_args(options);
    tm1Core.open_output_file();
    retVal = tm1Core.calculate_all_scores();
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

//
// TM1CoreInput methods
//
template<class TRNAString>
void TM1CoreInput<TRNAString>::init_from_args(TM1CSOptions &opts) {
    mMiRNAFasta = opts.mMiRNAFasta;
    mMRNAFasta = opts.mMRNAFasta;
}

template<class TRNAString>
int TM1CoreInput<TRNAString>::load_seq_from_file() {
    int retVal;

    // Read miRNA fasta file
    retVal = mMiRNASeqs.read_fasta(mMiRNAFasta);
    if (retVal != 0) {
        return 1;
    }

    // Read mRNA fasta file
    retVal = mMRNASeqs.read_fasta(mMRNAFasta);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

//
// TM1Core methods
//
template<class TRNAString, int SEEDLEN>
void TM1Core<TRNAString, SEEDLEN>::init_from_args(TM1CSOptions &opts) {
    mOutputAlign = opts.mOutputAlign;
    mOFileSite = opts.mOFileSite;
    mOFileScore = opts.mOFileScore;
}

template<class TRNAString, int SEEDLEN>
int TM1Core<TRNAString, SEEDLEN>::open_output_file() {
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

template<class TRNAString, int SEEDLEN>
int TM1Core<TRNAString, SEEDLEN>::calculate_all_scores() {
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
int TM1Core<TRNAString, SEEDLEN>::calculate_mirna_scores(unsigned pIdx) {
    int retVal;

    // Search seed sites
    if (mExecSearchSeedSites) {
        retVal = mSeedSites.find_seed_sites(mMiRNASeqs[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Seed site search failed." << std::endl;
            return 1;
        }
    }

    // Sort target sites
    if (mExecSortSites) {
        retVal = mSortedSites.generate_sorted_mrna_pos(mSeedSites, true);
        if (retVal != 0) {
            std::cerr << "ERROR: Sort target sites failed." << std::endl;
            return 1;
        }
    }

    // Get raw features
    if (mExecGetRawFeat) {
        retVal = mRawFeatures.add_features(mMiRNASeqs[pIdx], mMRNASeqs, mSeedSites, mSortedSites);
        if (retVal != 0) {
            std::cerr << "ERROR: Feature calculation failed." << std::endl;
            return 1;
        }
    }

    // Sort target sites again after creating site features
    if (mExecSortSites) {
        mSortedSites.clear_site_pos();
        retVal = mSortedSites.generate_sorted_mrna_pos(mSeedSites, false);
        if (retVal != 0) {
            std::cerr << "ERROR: Sort target sites failed." << std::endl;
            return 1;
        }
    }

    // Get mRNA features
    if (mExecGetMRNAFeat) {
        retVal = mMRNAFeatures.add_features(mSeedSites, mRawFeatures, mSortedSites);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate total ddG scores failed." << std::endl;
            return 1;
        }
    }

    // Calculate mRNA SVM scores
    if (mExecRNAScore) {
        retVal = mMRNAInput.classify(mMRNAFeatures.get_scaled_feature());
        if (retVal != 0) {
            std::cerr << "ERROR: mRNA SVM classification failed." << std::endl;
            return 1;
        }
    }


    // Summarize classified results
    if (mExecSumScores) {
        retVal = mScores.calc_scores(mMRNAFeatures.get_site_counts(), mMRNAInput.get_scores());
        if (retVal != 0) {
            std::cerr << "ERROR: Summarizing classified results failed." << std::endl;
            return 1;
        }
    }

    // Write site positions
    if (mOutputSitePos) {
        retVal = write_site_positions(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write site positions." << std::endl;
            return 1;
        }
    }

    // Write total scores
    if (mOutputScore) {
        retVal = write_scores(mMiRNAIds[pIdx]);
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
    mRawFeatures.clear_features();
    mSortedSites.clear_site_pos();
    mMRNAFeatures.clear_features();
    mMRNAInput.clear_scores();
    mScores.clear_scores();

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TM1Core<TRNAString, SEEDLEN>::write_site_positions(seqan::CharString const &pMiRNAId) {
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
        mOFile1 << toCString(mRawFeatures.get_seed_type(i)) << "\t";
        mOFile1 << 0;
        mOFile1 << std::endl;

    }

    return 0;

}

template<class TRNAString, int SEEDLEN>
int TM1Core<TRNAString, SEEDLEN>::write_scores(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &scores = mScores.get_scores();
    const seqan::String<int> &predictions = mScores.get_labels();
    const seqan::String<unsigned> &siteNum = mScores.get_site_num();
    const seqan::String<unsigned> &mRNAIdMap = mSortedSites.get_mrna_ids();
    typedef std::multimap<double, unsigned>::iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;
    TItMap itPos;
    std::multimap<double, unsigned> sortedMRNAByScore;

    for (unsigned i = 0; i < length(mRNAIdMap); ++i) {
        sortedMRNAByScore.insert(TPosPair(-1.0f * scores[i], i));
    }

    for (itPos = sortedMRNAByScore.begin(); itPos != sortedMRNAByScore.end(); ++itPos) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) mMRNAIds[mRNAIdMap[(*itPos).second]]) << "\t";
        mOFile2 << scores[(*itPos).second] << "\t";
        mOFile2 << siteNum[(*itPos).second] << "\t";
        mOFile2 << predictions[(*itPos).second];
        mOFile2 << std::endl;
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TM1Core<TRNAString, SEEDLEN>::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &mSeedTypes = mSeedSites.get_seed_types();
    const TM1Alignment<TRNAString> &alignment = mRawFeatures.get_alignment();

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
        alignment.write_alignment(i);
        std::cout << "  miRNA:                " << toCString(pMiRNAId) << std::endl;
        std::cout << "  mRNA:                 " << toCString((seqan::CharString) mMRNAIds[mRNAPos[i]]) << std::endl;
        std::cout << "  seed type:            " << toCString(seedType) << std::endl;
        std::cout << "  seed type:            " << toCString(mRawFeatures.get_seed_type(i));
        std::cout << std::endl;
        std::cout << "  position(seed start): " << seedStart << std::endl;
        std::cout << "  position(seed end):   " << seedEnd << std::endl;
        std::cout << std::endl << std::endl;

        ++count;

    }

    return 0;
}

// Explicit template instantiation
template
class TM1CoreInput<TRNATYPE>;

template
class TM1Core<TRNATYPE>;

} // namespace tm1p
