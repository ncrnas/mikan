#include <mikan/lib/targetscan5_cs/include/ts5_core.hpp>          // TS5CoreInput, TS5Core
#include <mikan/lib/targetscan5_cs/include/ts5_inst_template.hpp> // TRNATYPE
#include <mikan/lib/targetscan5_cs/include/ts5_option.hpp>        // TS5CSOptions
#include <seqan/arg_parse.h>

int main(int argc, char const ** argv)
{
    int retVal;

    // Parse the command line
    ts5cs::TS5CSOptions options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
    {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    ts5cs::TS5CoreInput<ts5cs::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0)
    {
        return 1;
    }

    // Create index
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TRNASet const& mMRNASeqs = coreInput.get_mrna_seqs();
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TCharSet const& mMiRNAIds = coreInput.get_mirna_ids();
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TRNASet const& mMiRNASeqs = coreInput.get_mirna_seqs();
    ts5cs::TS5Core<ts5cs::TRNATYPE>::TCharSet const& mMRNAIds = coreInput.get_mrna_ids();
    ts5cs::TS5Core<ts5cs::TRNATYPE> ts5Core(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    ts5Core.init_from_args(options);
    ts5Core.open_output_file();
    retVal = ts5Core.calculate_all_scores();
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
