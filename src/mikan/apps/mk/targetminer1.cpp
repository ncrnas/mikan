#include <mikan/lib/targetminer1/include/tm1_core.hpp>          // TM1CoreInput, TM1Core

int main(int argc, char const ** argv)
{
    int retVal;

    // Parse the command line
    tm1p::TM1CSOptions options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
    {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    tm1p::TM1CoreInput<tm1p::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0)
    {
        return 1;
    }

    // Create index
    tm1p::TM1Core<tm1p::TRNATYPE>::TRNASet const& mMRNASeqs = coreInput.get_mrna_seqs();
    tm1p::TM1Core<tm1p::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    tm1p::TM1Core<tm1p::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    tm1p::TM1Core<tm1p::TRNATYPE>::TCharSet const& mMiRNAIds = coreInput.get_mirna_ids();
    tm1p::TM1Core<tm1p::TRNATYPE>::TRNASet const& mMiRNASeqs = coreInput.get_mirna_seqs();
    tm1p::TM1Core<tm1p::TRNATYPE>::TCharSet const& mMRNAIds = coreInput.get_mrna_ids();
    tm1p::TM1Core<tm1p::TRNATYPE> tm1Core(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    tm1Core.init_from_args(options);
    tm1Core.open_output_file();
    retVal = tm1Core.calculate_all_scores();
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
