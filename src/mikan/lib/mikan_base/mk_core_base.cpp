#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include "mk_core_base.hpp"      // MKCoreBase

namespace mikan {

//
// MKCoreBase methods
//
void MKCoreBase::init_from_args(mikan::MKOptions const &opts) {
    mOutputAlign = opts.mOutputAlign;
    mOFileSite = opts.mOFileSite;
    mOFileRNA = opts.mOFileTotal;
}

int MKCoreBase::open_output_file() {
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

void MKCoreBase::close_output_file() {
    // Close output file 1
    if (mOFile1.is_open()) {
        mOFile1.close();
    }

    // Close output file 2
    if (mOFile2.is_open()) {
        mOFile2.close();
    }
}

int MKCoreBase::calculate_all_scores() {
    int retVal;

    for (unsigned i = 0; i < length(mMiRNASeqs); ++i) {

#if SEQAN_ENABLE_DEBUG
        clock_t startTime = clock();
#endif

        retVal = calculate_mirna_scores(i);
        if (retVal != 0) {
            std::cerr << "ERROR: Score calculation failed for ";
            std::cerr << toCString((mikan::TCharStr) mMiRNAIds[i]) << "." << std::endl;
            return 1;
        }

#if SEQAN_ENABLE_DEBUG
        std::cout << toCString((mikan::TCharStr) mMiRNAIds[i]) << ": ";
        std::cout << double(clock() - startTime) / (double) CLOCKS_PER_SEC << " seconds." << std::endl;
#endif

    }

    return 0;
}

int MKCoreBase::calculate_mirna_scores(unsigned pIdx) {
    int retVal;

    retVal = find_seed_sites(pIdx);
    if (retVal != 0) {
        std::cerr << "ERROR: Finding seed sites failed." << std::endl;
        return 1;
    }

    retVal = calc_site_scores(pIdx);
    if (retVal != 0) {
        std::cerr << "ERROR: Calculating seed scores failed." << std::endl;
        return 1;
    }

    retVal = calc_rna_scores(pIdx);
    if (retVal != 0) {
        std::cerr << "ERROR: Calculating RNA scores failed." << std::endl;
        return 1;
    }

    retVal = output_results(pIdx);
    if (retVal != 0) {
        std::cerr << "ERROR: Output results failed." << std::endl;
        return 1;
    }

    clear_all();

    return 0;
}

int MKCoreBase::write_site_score(
        mikan::TCharStr const &pMiRNAId,
        mikan::MKSeedSites &pSeedSites,
        mikan::MKRMAWithSites &pRNAWithSites) {

    seqan::StringSet<seqan::String<unsigned> > &rnaSitePosMap = pRNAWithSites.get_rna_site_pos_map();
    mikan::TMRNAPosSet &uniqRNAPosSet = pRNAWithSites.get_uniq_mrna_pos_set();

    int retVal = 0;
    int count = 0;
    for (unsigned i = 0; i < length(pRNAWithSites.mEffectiveRNAs); i++) {
        if (!pRNAWithSites.mEffectiveRNAs[i]) {
            continue;
        }

        for (unsigned j = 0; j < length(rnaSitePosMap[i]); ++j) {
            if (!pSeedSites.mEffectiveSites[rnaSitePosMap[i][j]]) {
                continue;
            }
            prepare_site_output(pMiRNAId, uniqRNAPosSet[i], rnaSitePosMap[i][j]);
        }

        ++count;
    }

    return retVal;
}

void MKCoreBase::write_site_score_tab(
        std::string &pMiRNAName,
        std::string &pMRNAName,
        unsigned pStartPos,
        unsigned pEndPos,
        std::string &pSeedType,
        std::string &pScore1Name,
        std::string &pScore1,
        std::string &pScore2Name,
        std::string &pScore2) {

    if (mPrintSiteHeader) {
        mOFile1 << "# miRNA name, ";
        mOFile1 << "mRNA name, ";
        mOFile1 << "start (1-base), ";
        mOFile1 << "end (1-base), ";
        mOFile1 << "seed type, ";
        mOFile1 << "score 1 (" << pScore1Name << "), ";
        mOFile1 << "score 2 (" << pScore2Name << ")";
        mOFile1 << std::endl;
        mPrintSiteHeader = false;
    }

    mOFile1 << pMiRNAName << "\t";
    mOFile1 << pMRNAName << "\t";
    mOFile1 << pStartPos << "\t";
    mOFile1 << pEndPos << "\t";
    mOFile1 << pSeedType << "\t";
    mOFile1 << pScore1 << "\t";
    mOFile1 << pScore2;
    mOFile1 << std::endl;

}

int MKCoreBase::write_rna_score(mikan::TCharStr const &pMiRNAId) {
    int retVal = 0;

    prepare_rna_output(pMiRNAId);

    return retVal;
}

void MKCoreBase::write_rna_score_tab(
        std::string &pMiRNAName,
        std::string &pMRNAName,
        unsigned pSiteNum,
        std::string &pScore1Name,
        std::string &pScore1,
        std::string &pScore2Name,
        std::string &pScore2) {

    if (mPrintRNAheader) {
        mOFile2 << "# miRNA name, ";
        mOFile2 << "mRNA name, ";
        mOFile2 << "number of sites, ";
        mOFile2 << "score 1 (" << pScore1Name << "), ";
        mOFile2 << "score 2 (" << pScore2Name << ")";
        mOFile2 << std::endl;
        mPrintRNAheader = false;
    }


    mOFile2 << pMiRNAName << "\t";
    mOFile2 << pMRNAName << "\t";
    mOFile2 << pSiteNum << "\t";
    mOFile2 << pScore1 << "\t";
    mOFile2 << pScore2;
    mOFile2 << std::endl;

}

} // namespace mr3as
