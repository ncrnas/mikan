#include <cmath>                  // roundf
#include <iostream>
//#define SEQAN_ENABLE_DEBUG 1
#if SEQAN_ENABLE_DEBUG
#include <ctime>                  // clock_t, clock, CLOCKS_PER_SEC
#endif

#include "mk_core.hpp"            // MKCoreBase

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

int MKCoreBase::calculate_all_scores() {
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

    retVal = ensemble_site_scores(pIdx);
    if (retVal != 0) {
        std::cerr << "ERROR: Ensemble site scores failed." << std::endl;
        return 1;
    }

    retVal = calc_rna_scores(pIdx);
    if (retVal != 0) {
        std::cerr << "ERROR: Calculating RNA scores failed." << std::endl;
        return 1;
    }

    retVal = ensemble_rna_scores(pIdx);
    if (retVal != 0) {
        std::cerr << "ERROR: Ensemble RNA scores failed." << std::endl;
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

} // namespace mr3as
