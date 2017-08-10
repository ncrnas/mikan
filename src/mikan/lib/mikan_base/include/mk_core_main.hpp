#ifndef MK_CORE_MAIN_HPP_
#define MK_CORE_MAIN_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_input.hpp"           // MKInput

namespace mikan {

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


} // namespace mikan

#endif /* MK_CORE_MAIN_HPP_ */

