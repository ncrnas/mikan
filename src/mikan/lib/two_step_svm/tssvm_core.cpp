#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                    // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_inst_template.hpp"     // TRNATYPE
#include "tssvm_option.hpp"         // TSSVMOptions
#include "tssvm_seed_site.hpp"      // TSSVMSeedSites, TSSVMSeedSiteOverlap
#include "tssvm_align.hpp"          // TSAlign
#include "tssvm_site_feature.hpp"   // TSSVMRawFeatures
#include "tssvm_site_svm.hpp"       // TSSVMSiteModel, TSSVMSiteInputVector
#include "tssvm_mrna_feature.hpp"   // TSSVMRNARawFeatures
#include "tssvm_mrna_svm.hpp"       // TSSVMRNAInputVector
#include "tssvm_core.hpp"           // TSSVMCore
#include "mk_input.hpp"             // MKInput

using namespace mikan;

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
    mikan::MKInput<TRNATYPE> coreInput;
    coreInput.set_file_names(options.mMiRNAFasta, options.mMRNAFasta);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0) {
        return 1;
    }

    // Create index
    tssvm::TSSVMCore<TRNATYPE>::TRNASet const &mMRNASeqs = coreInput.get_mrna_seqs();
    tssvm::TSSVMCore<TRNATYPE>::TIndexQGram index(mMRNASeqs);
    tssvm::TSSVMCore<TRNATYPE>::TFinder finder(index);

    // Prepare models
    tssvm::TSSVMCore<TRNATYPE>::TCharSet const &mMiRNAIds = coreInput.get_mirna_ids();
    tssvm::TSSVMCore<TRNATYPE>::TRNASet const &mMiRNASeqs = coreInput.get_mirna_seqs();
    tssvm::TSSVMCore<TRNATYPE>::TCharSet const &mMRNAIds = coreInput.get_mrna_ids();
    tssvm::TSSVMCore<TRNATYPE> tssvmCore(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    tssvmCore.init_from_args(options);
    retVal = tssvmCore.init_site_svm();
    if (retVal != 0) {
        return 1;
    }

    // Calculate scores for all miRNAs
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
template<class TRNAString, int SEEDLEN>
void TSSVMCore<TRNAString, SEEDLEN>::init_from_args(TSSVMOptions &opts) {
    mOFileTargetSite = opts.mOFileSite;
    mOFileMRNA = opts.mOFileTotal;
    mOutputAlign = opts.mOutputAlign;
}

template<class TRNAString, int SEEDLEN>
int TSSVMCore<TRNAString, SEEDLEN>::init_site_svm() {
    return mSiteModel.init_model();
}

template<class TRNAString, int SEEDLEN>
int TSSVMCore<TRNAString, SEEDLEN>::open_output_file() {
    // Open output file 1
    mOFile1.open(toCString(mOFileTargetSite), std::ofstream::out);
    if (!mOFile1.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileTargetSite) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    // Open output file 2
    mOFile2.open(toCString(mOFileMRNA), std::ofstream::out);
    if (!mOFile2.good()) {
        std::cerr << "ERROR: Could not open output file " << toCString(mOFileMRNA) << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TSSVMCore<TRNAString, SEEDLEN>::calculate_all_scores() {
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
int TSSVMCore<TRNAString, SEEDLEN>::calculate_mirna_scores(unsigned pIdx) {
    int retVal;
    TRNAString miRNASeq = mMiRNASeqs[pIdx];

    // Search seed sites
    if (mExecSearchSeedSites) {
        retVal = mSeedSites.find_seed_sites(miRNASeq);
        if (retVal != 0) {
            std::cerr << "ERROR: Seed site search failed." << std::endl;
            return 1;
        }
    }

    // Filter overlapped sites
    if (mExecFilterOverlap) {
        retVal = mOverlappedSites.filter_overlapped_sites(mSeedSites, (unsigned) length(mMRNASeqs));
        if (retVal != 0) {
            std::cerr << "ERROR: Check overlapped sites failed." << std::endl;
            return 1;
        }
    }

    // Align sequences
    if (mExecAlignSeq) {
        retVal = mAlignSeqs.align_seq(mSeedSites, miRNASeq, mMRNASeqs);
        if (retVal != 0) {
            std::cerr << "ERROR: Align sequences failed." << std::endl;
            return 1;
        }
    }

    // Generate site features
    if (mExecSiteFeat) {
        retVal = mSiteFeatures.add_features(mSeedSites, mAlignSeqs, miRNASeq, mMRNASeqs);
        if (retVal != 0) {
            std::cerr << "ERROR: Site feature calculation failed." << std::endl;
            return 1;
        }
    }

    // Calculate site SVM scores
    if (mExecSiteScore) {
        retVal = mSiteInput.classify(mSiteFeatures);
        if (retVal != 0) {
            std::cerr << "ERROR: Site SVM classification failed." << std::endl;
            return 1;
        }
    }

    // Generate RNA features
    if (mExecRNAFeat) {
        retVal = mRnaFeatures.add_features(mSeedSites, mMRNASeqs, mOverlappedSites, mSiteInput);
        if (retVal != 0) {
            std::cerr << "ERROR: RNA feature calculation failed." << std::endl;
            return 1;
        }
    }

    // Calculate RNA SVM scores
    if (mExecRNAScore) {
        retVal = mRnaInput.classify(mRnaFeatures);
        if (retVal != 0) {
            std::cerr << "ERROR: RNA SVM classification failed." << std::endl;
            return 1;
        }
    }

    // Write TargetSite scores
    if (mOutputSiteScore) {
        retVal = write_ts_scores(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            std::cerr << "ERROR: Could not write target-site scores." << std::endl;
            return 1;
        }
    }

    // Write mRNA scores
    if (mOutputRNAScore) {
        retVal = write_mrna_scores(mMiRNAIds[pIdx]);
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

    mSeedSites.clear_pos();
    mOverlappedSites.clear_site_pos();
    mAlignSeqs.clear_alignments();
    mSiteFeatures.clear_features();
    mSiteInput.clear_scores();
    mRnaFeatures.clear_features();
    mRnaInput.clear_scores();

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TSSVMCore<TRNAString, SEEDLEN>::write_ts_scores(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();
    seqan::CharString seedType;
    int seedStart;
    std::set<unsigned>::iterator itSet;
    std::set<unsigned> &rnaPosSet = mOverlappedSites.get_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &sortedMRNAPos = mOverlappedSites.get_sorted_mrna_pos();
    const seqan::String<float> &scors = mSiteInput.get_scores();

    for (itSet = rnaPosSet.begin(); itSet != rnaPosSet.end(); ++itSet) {
        for (unsigned i = 0; i < length(sortedMRNAPos[*itSet]); ++i) {
            if (!mSeedSites.mEffectiveSites[sortedMRNAPos[*itSet][i]]) {
                continue;
            }

            seedType = seedTypes[sortedMRNAPos[*itSet][i]];
            seedStart = sitePos[sortedMRNAPos[*itSet][i]];
            if (seedType == "7mer-A1") {
                seedStart += 1;
            }

            mOFile1 << toCString(pMiRNAId) << "\t";
            mOFile1 << toCString((seqan::CharString) mMRNAIds[mRNAPos[sortedMRNAPos[*itSet][i]]]) << "\t";
            mOFile1 << seedStart + 1 << "\t";
            mOFile1 << seedStart + 7 << "\t";
            mOFile1 << toCString((seqan::CharString) seedTypes[sortedMRNAPos[*itSet][i]]) << "\t";
            mOFile1 << scors[sortedMRNAPos[*itSet][i]] << "\t";
            mOFile1 << std::endl;
        }
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TSSVMCore<TRNAString, SEEDLEN>::write_mrna_scores(seqan::CharString const &pMiRNAId) {
    typedef std::multimap<float, unsigned>::reverse_iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;
    std::set<unsigned>::iterator itSet;
    TItMap itPos;
    std::set<unsigned> &rnaPosSet = mOverlappedSites.get_mrna_pos_set();
    const seqan::String<float> &scors = mRnaInput.get_scores();
    std::multimap<float, unsigned> sortedMRNAByScore;
    seqan::String<unsigned> &siteCount = mRnaFeatures.get_site_count();

    for (itSet = rnaPosSet.begin(); itSet != rnaPosSet.end(); ++itSet) {

        if (!mRnaFeatures.mEffectiveRNAs[*itSet]) {
            continue;
        }

        sortedMRNAByScore.insert(TPosPair((float) scors[*itSet], *itSet));

    }

    for (itPos = sortedMRNAByScore.rbegin(); itPos != sortedMRNAByScore.rend(); ++itPos) {
        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) mMRNAIds[(*itPos).second]) << "\t";
        mOFile2 << scors[(*itPos).second] << "\t";
        mOFile2 << siteCount[(*itPos).second] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

template<class TRNAString, int SEEDLEN>
int TSSVMCore<TRNAString, SEEDLEN>::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();
    seqan::CharString seedType;
    int seedStart;
    std::set<unsigned>::iterator itSet;
    std::set<unsigned> &rnaPosSet = mOverlappedSites.get_mrna_pos_set();
    seqan::StringSet<seqan::String<unsigned> > &sortedMRNAPos = mOverlappedSites.get_sorted_mrna_pos();
    const seqan::String<float> &scors = mSiteInput.get_scores();
    int count = 0;

    for (itSet = rnaPosSet.begin(); itSet != rnaPosSet.end(); ++itSet) {
        for (unsigned i = 0; i < seqan::length(sortedMRNAPos[*itSet]); ++i) {
            if (!mSeedSites.mEffectiveSites[sortedMRNAPos[*itSet][i]]) {
                continue;
            }

            seedType = seedTypes[sortedMRNAPos[*itSet][i]];
            seedStart = sitePos[sortedMRNAPos[*itSet][i]];
            if (seedType == "7mer-A1") {
                seedStart += 1;
            }

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mAlignSeqs.write_alignment(sortedMRNAPos[*itSet][i]);
            std::cout << "  miRNA:                " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:                 ";
            std::cout << toCString((seqan::CharString) mMRNAIds[mRNAPos[sortedMRNAPos[*itSet][i]]]) << std::endl;
            std::cout << "  seed type:            " << toCString(seedType) << std::endl;
            std::cout << "  position(seed start): " << seedStart + 1 << std::endl;
            std::cout << "  site level score:     " << scors[sortedMRNAPos[*itSet][i]];
            std::cout << std::endl << std::endl;

            ++count;
        }
    }

    return 0;
}

// Explicit template instantiation

template
class TSSVMCore<TRNATYPE>;

} // namespace tssvm
