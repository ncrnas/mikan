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
int MR3Core::find_seed_sites(unsigned pIdx) {
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

int MR3Core::calc_site_scores(unsigned pIdx) {
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

    return 0;

}

int MR3Core::ensemble_site_scores(unsigned) {

    return 0;

}

int MR3Core::calc_rna_scores(unsigned) {
    int retVal;

    if (mCalcRNAScore) {
        retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
        if (retVal != 0) {
            return 1;
        }
    }

    return 0;

}

int MR3Core::ensemble_rna_scores(unsigned) {

    return 0;

}

int MR3Core::output_results(unsigned pIdx) {
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

void MR3Core::clear_all() {
    mSeedSeqs.clear_seeds();
    mSeedSites.clear_pos();
    mSiteScores.clear_scores();
    mRNAWithSites.clear_maps();
    mRNAScores.clear_scores();
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
