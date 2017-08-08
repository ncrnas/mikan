#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                    // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"             // MKInput
#include "tssvm_option.hpp"         // TSSVMOptions
#include "tssvm_seed_site.hpp"      // TSSVMSeedSites, TSSVMSiteFilter
#include "tssvm_align.hpp"          // TSAlign
#include "tssvm_site_feature.hpp"   // TSSVMRawFeatures
#include "tssvm_site_svm.hpp"       // TSSVMSiteModel, TSSVMSiteInputVector
#include "tssvm_mrna_feature.hpp"   // TSSVMRNARawFeatures
#include "tssvm_mrna_svm.hpp"       // TSSVMRNAInputVector
#include "tssvm_core.hpp"           // TSSVMCore

namespace tssvm {

int TSSVMCoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line.
    tssvm::TSSVMOptions options;
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
    tssvm::TSSVMCore tssvmCore(options, mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    tssvmCore.open_output_file();
    retVal = tssvmCore.calculate_all_scores();
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

//
// TSSVMCore methods
//
void TSSVMCore::init_from_args(mikan::MKOptions const &opts) {
    mOFileSite = opts.mOFileSite;
    mOFileRNA = opts.mOFileTotal;
    mOutputAlign = opts.mOutputAlign;

    resize(mSeedTypeDef, 1);
    mSeedTypeDef[0] = "";
}

int TSSVMCore::open_output_file() {
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

int TSSVMCore::calculate_all_scores() {
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

int TSSVMCore::calculate_mirna_scores(unsigned pIdx) {
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

    // Filter overlapped sites
    mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
    if (mExecFilterOverlap) {
        retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            std::cerr << "ERROR: Check overlapped sites failed." << std::endl;
            return 1;
        }
    }

    // Calculate site SVM scores
    if (mExecSiteScore) {
        retVal = mSiteScores.calc_scores(miRNASeq, mMRNASeqs, mSeedSites, mRNAWithSites);
        if (retVal != 0) {
            std::cerr << "ERROR: Calculate site SVM scores failed." << std::endl;
            return 1;
        }
    }

    // Calculate RNA SVM scores
    if (mExecRNAScore) {
        retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            std::cerr << "ERROR: RNA SVM classification failed." << std::endl;
            return 1;
        }
    }

    // Write TargetSite scores
    if (mOutputSiteScore) {
        retVal = write_site_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write target-site scores." << std::endl;
            return 1;
        }
    }

    // Write mRNA scores
    if (mOutputRNAScore) {
        retVal = write_rna_score(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write mRNA scores." << std::endl;
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
    mRNAWithSites.clear_maps();
    mSiteScores.clear_scores();
    mRNAScores.clear_scores();

    return 0;
}

int TSSVMCore::write_site_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<float> &scores = mSiteScores.get_scores();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();

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
            if (seedType == "7mer-A1") {
                seedStart += 1;
            }

            float score = scores[rnaSitePosMap[i][j]];
            score = roundf(score * 10000.0f) / 10000.0f;
            mOFile1 << toCString(pMiRNAId) << "\t";
            mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[rnaSitePosMap[i][j]]]) << "\t";
            mOFile1 << seedStart + 1 << "\t";
            mOFile1 << seedStart + 7 << "\t";
            mOFile1 << toCString((seqan::CharString) seedTypes[rnaSitePosMap[i][j]]) << "\t";
            mOFile1 << score << "\t";
            mOFile1 << std::endl;
        }

    }

    return 0;

}

int TSSVMCore::write_rna_score(seqan::CharString const &pMiRNAId) {
    typedef std::multimap<float, unsigned>::reverse_iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;

    const seqan::String<float> &scores = mRNAScores.get_scores();
    const seqan::String<int> &siteCount = mRNAScores.get_site_num();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

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
        mOFile2 << toCString((seqan::CharString) mMRNAIds[uniqRNAPosSet[(*itPos).second]]) << "\t";
        mOFile2 << score << "\t";
        mOFile2 << siteCount[(*itPos).second] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

int TSSVMCore::write_alignment(seqan::CharString const &pMiRNAId) {

    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();
    const seqan::String<float> &scores = mSiteScores.get_scores();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();

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
            if (seedType == "7mer-A1") {
                seedStart += 1;
            }

            std::cout << "### " << (i + j) + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.write_alignment(rnaSitePosMap[i][j]);
            std::cout << "  miRNA:                " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:                 ";
            std::cout << toCString((seqan::CharString) mMRNAIds[mRNAPos[rnaSitePosMap[i][j]]]) << std::endl;
            std::cout << "  seed type:            " << toCString(seedType) << std::endl;
            std::cout << "  position(seed start): " << seedStart + 1 << std::endl;
            std::cout << "  site level score:     " << scores[rnaSitePosMap[i][j]];
            std::cout << std::endl << std::endl;


        }
    }

    return 0;

}

} // namespace tssvm
