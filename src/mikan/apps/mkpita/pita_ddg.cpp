#include <mikan/lib/pita_ddg/include/pita_core.hpp>          // PITACoreInput, PITACore

int main(int argc, char const ** argv)
{
    int retVal;

    // Parse the command line.
    ptddg::PITAOptions options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
    {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    ptddg::PITACoreInput<ptddg::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0)
    {
        return 1;
    }

    // Create index
    ptddg::PITACore<ptddg::TRNATYPE>::TRNASet const& mMRNASeqs = coreInput.get_mrna_seqs();
    ptddg::PITACore<ptddg::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    ptddg::PITACore<ptddg::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    ptddg::PITACore<ptddg::TRNATYPE>::TCharSet const& mMiRNAIds = coreInput.get_mirna_ids();
    ptddg::PITACore<ptddg::TRNATYPE>::TRNASet const& mMiRNASeqs = coreInput.get_mirna_seqs();
    ptddg::PITACore<ptddg::TRNATYPE>::TCharSet const& mMRNAIds = coreInput.get_mrna_ids();

    ptddg::PITACore<ptddg::TRNATYPE> pitaCore(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    pitaCore.init_from_args(options);
    pitaCore.open_output_file();
    retVal = pitaCore.calculate_all_scores();
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
