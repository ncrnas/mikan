#ifndef MIKAN_TEST_MKENS_HPP_
#define MIKAN_TEST_MKENS_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mk_option.hpp"

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mk_main.hpp"
#include "mke_core.hpp"
#include "mk_input.hpp"
#include "mk_seed_site.hpp"

typedef TestSite<mkens::MKEOptions, mikan::MKSeedSeqs, mikan::MKSeedSites> TestSiteMKENS;

template<class TOptions, class TCore>
int MKCoreMain(int argc, char const **argv) {
    int retVal;

    // Parse the command line.
    TOptions options;
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

    TCore mkCore(options, mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    mkCore.open_output_file();
    retVal = mkCore.calculate_all_scores();

    return retVal;
}

#endif //MIKAN_TEST_MKENS_HPP_
