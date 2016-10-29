#include <mr3_core.hpp>          // MR3CoreInput, MR3Core

int main(int argc, char const ** argv)
{
    int retVal;

    // Parse the command line.
    mr3as::MR3Options options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
    {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    mr3as::MR3CoreInput<mr3as::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0)
    {
        return 1;
    }

    // Create index
    mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet const& mMRNASeqs = coreInput.get_mrna_seqs();
    mr3as::MR3Core<mr3as::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    mr3as::MR3Core<mr3as::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet const& mMiRNAIds = coreInput.get_mirna_ids();
    mr3as::MR3Core<mr3as::TRNATYPE>::TRNASet const& mMiRNASeqs = coreInput.get_mirna_seqs();
    mr3as::MR3Core<mr3as::TRNATYPE>::TCharSet const& mMRNAIds = coreInput.get_mrna_ids();

    mr3as::MR3Core<mr3as::TRNATYPE> pitaCore(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    pitaCore.init_from_args(options);
    pitaCore.open_output_file();
    retVal = pitaCore.calculate_all_scores();
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
