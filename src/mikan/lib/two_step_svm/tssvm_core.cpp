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
int TSSVMCore::find_seed_sites(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    if (mFindSeedSites) {
        retVal = mSeedSeqs.create_seed_seqs(miRNASeq);
        if (retVal != 0) {
            return 1;
        }

        retVal = mSeedSites.find_seed_sites(mSeedSeqs);
        if (retVal != 0) {
            return 1;
        }

    }

    if (mFilterSites) {
        mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
        retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;
}

int TSSVMCore::calc_site_scores(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    if (mCalcSiteScore) {
        retVal = mSiteScores.calc_scores(miRNASeq, mMRNASeqs, mSeedSites, mRNAWithSites);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;

}

int TSSVMCore::ensemble_site_scores(unsigned) {

    return 0;

}

int TSSVMCore::calc_rna_scores(unsigned) {
    int retVal;

    if (mCalcRNAScore) {
        mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
        retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;

}

int TSSVMCore::ensemble_rna_scores(unsigned) {

    return 0;

}

int TSSVMCore::output_results(unsigned pIdx) {
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

    // Write alignments
    if (mOutputAlign) {
        retVal = write_alignment(mMiRNAIds[pIdx]);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;

}

void TSSVMCore::clear_all() {
    mSeedSeqs.clear_seeds();
    mSeedSites.clear_pos();
    mRNAWithSites.clear_maps();
    mSiteScores.clear_scores();
    mRNAScores.clear_scores();
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
