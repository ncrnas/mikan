#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include <seqan/arg_parse.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_input.hpp"           // MKInput
#include "pita_option.hpp"        // PITAOptions
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_score.hpp"    // PITAMFEScores
#include "pita_site_filter.hpp"   // PITASiteFilter
#include "pita_core.hpp"          // PITACore

namespace ptddg {

int PITACoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line.
    ptddg::PITAOptions options;
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

    ptddg::PITACore pitaCore(options, mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    pitaCore.open_output_file();
    retVal = pitaCore.calculate_all_scores();
    if (retVal != 0) {
        return 1;
    }

    return 0;
}

//
// PITACore methods
//
int PITACore::find_seed_sites(unsigned pIdx) {
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

int PITACore::calc_site_scores(unsigned pIdx) {
    int retVal;
    mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

    if (mCalcSiteScore) {
        retVal = mSiteScores.calc_scores(miRNASeq, mMRNASeqs, mSeedSites, mRNAWithSites);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mFilterSiteScores) {
        mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
        retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            return 1;
        }
    }

    if (mSelectTopSites) {
        retVal = mTopNSites.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;

}

int PITACore::ensemble_site_scores(unsigned) {

    return 0;

}

int PITACore::calc_rna_scores(unsigned) {
    int retVal;

    if (mCalcRNAScore) {
        retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;

}

int PITACore::ensemble_rna_scores(unsigned) {

    return 0;

}

int PITACore::output_results(unsigned pIdx) {
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

void PITACore::clear_all() {
    mSeedSeqs.clear_seeds();
    mSeedSites.clear_pos();
    mRNAWithSites.clear_maps();
    mSiteScores.clear_scores();
    mRNAScores.clear_scores();
}

int PITACore::write_site_score(seqan::CharString const &pMiRNAId) {
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
            score = roundf(score * 100.0f) / 100.0f;

            mOFile1 << toCString(pMiRNAId) << "\t";
            mOFile1 << toCString((seqan::CharString) (mMRNAIds[uniqRNAPosSet[i]])) << "\t";
            mOFile1 << seedStart + 1 << "\t";
            mOFile1 << seedStart + 1 + mikan::SEEDLEN << "\t";
            mOFile1 << toCString((seqan::CharString) (seedTypes[rnaSitePosMap[i][j]])) << "\t";
            mOFile1 << score << "\t";
            mOFile1 << std::endl;

            ++count;
        }

    }

    return 0;

}

int PITACore::write_rna_score(seqan::CharString const &pMiRNAId) {
    const seqan::String<float> &totalScores = mRNAScores.get_scores();
    const seqan::String<int> &mRNAPos = mRNAScores.get_mrna_pos();
    const seqan::String<int> &siteNum = mRNAScores.get_site_num();
    float score;

    for (unsigned i = 0; i < length(mRNAPos); ++i) {
        score = totalScores[i];
        score = roundf(score * 100.0f) / 100.0f;

        mOFile2 << toCString(pMiRNAId) << "\t";
        mOFile2 << toCString((seqan::CharString) (mMRNAIds[mRNAPos[i]])) << "\t";
        mOFile2 << score << "\t";
        mOFile2 << siteNum[i] << "\t";
        mOFile2 << std::endl;
    }

    return 0;
}

int PITACore::write_alignment(seqan::CharString const &pMiRNAId) {
    const seqan::String<unsigned> &mRNAPos = mSeedSites.get_mrna_pos();
    const seqan::String<unsigned> &sitePos = mSeedSites.get_site_pos();
    const seqan::StringSet<seqan::CharString> &seedTypes = mSeedSites.get_seed_types();

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = mRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = mRNAWithSites.get_uniq_mrna_pos_set();

    seqan::CharString seedType;
    float dGduplex;
    float dG5;
    float dG3;
    float dGopen;
    float dG0;
    float dG1;
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
            score = roundf(score * 100.0f) / 100.0f;

            dGduplex = (float) mSiteScores.get_dgall(rnaSitePosMap[i][j]);
            dGduplex = roundf(dGduplex * 100.0f) / 100.0f;
            dG5 = (float) mSiteScores.get_dg5(rnaSitePosMap[i][j]);
            dG5 = roundf(dG5 * 100.0f) / 100.0f;
            dG3 = (float) mSiteScores.get_dg3(rnaSitePosMap[i][j]);
            dG3 = roundf(dG3 * 100.0f) / 100.0f;
            dGopen =
                    (float) mSiteScores.get_dg0(rnaSitePosMap[i][j]) - (float) mSiteScores.get_dg1(rnaSitePosMap[i][j]);
            dGopen = roundf(dGopen * 100.0f) / 100.0f;
            dG0 = (float) mSiteScores.get_dg0(rnaSitePosMap[i][j]);
            dG0 = roundf(dG0 * 100.0f) / 100.0f;
            dG1 = (float) mSiteScores.get_dg1(rnaSitePosMap[i][j]);
            dG1 = roundf(dG1 * 100.0f) / 100.0f;

            std::cout << "### " << count + 1 << ": " << toCString(pMiRNAId) << " ###" << std::endl;
            mSiteScores.print_alignment(rnaSitePosMap[i][j]);
            std::cout << "  miRNA:               " << toCString(pMiRNAId) << std::endl;
            std::cout << "  mRNA:                "
                      << toCString((seqan::CharString) (mMRNAIds[mRNAPos[uniqRNAPosSet[i]]]));
            std::cout << std::endl;
            std::cout << "  seed type:           " << toCString((seqan::CharString) (seedTypes[rnaSitePosMap[i][j]]))
                      << std::endl;
            std::cout << "  position(start):     " << seedStart + 1 << std::endl;
            std::cout << "  position(end):       " << seedStart + 1 + mikan::SEEDLEN << std::endl;
            std::cout << "  ddG:                 " << score << std::endl;
            std::cout << "  dGduplex(dG5 + dG3): " << dGduplex << std::endl;
            std::cout << "  dG5:                 " << dG5 << std::endl;
            std::cout << "  dG3:                 " << dG3 << std::endl;
            std::cout << "  dGopen(dG0 - dG1):   " << dGopen << std::endl;
            std::cout << "  dG0:                 " << dG0 << std::endl;
            std::cout << "  dG1:                 " << dG1 << std::endl;
            std::cout << std::endl;
        }

    }

    return 0;

}

} // namespace ptddg

