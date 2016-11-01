#include <tssvm_core.hpp>           // TSSVMCoreInput, TSSVMCore
#include <tssvm_inst_template.hpp>  // TRNATYPE

int main(int argc, char const ** argv)
{
    int retVal;

    // Parse the command line.
    tssvm::TSSVMOptions options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
    {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    tssvm::TSSVMCoreInput<tssvm::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0)
    {
        return 1;
    }

    // Create index
    tssvm::TSSVMCore<tssvm::TRNATYPE>::TRNASet const& mMRNASeqs = coreInput.get_mrna_seqs();
    tssvm::TSSVMCore<tssvm::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    tssvm::TSSVMCore<tssvm::TRNATYPE>::TFinder finder(index);

    // Prepare models
    tssvm::TSSVMCore<tssvm::TRNATYPE>::TCharSet const& mMiRNAIds = coreInput.get_mirna_ids();
    tssvm::TSSVMCore<tssvm::TRNATYPE>::TRNASet const& mMiRNASeqs = coreInput.get_mirna_seqs();
    tssvm::TSSVMCore<tssvm::TRNATYPE>::TCharSet const& mMRNAIds = coreInput.get_mrna_ids();
    tssvm::TSSVMCore<tssvm::TRNATYPE> tssvmCore(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder);
    tssvmCore.init_from_args(options);
    retVal = tssvmCore.init_site_svm();
    if (retVal != 0)
    {
        return 1;
    }

    // Calculate scores for all miRNAs
    tssvmCore.open_output_file();
    retVal = tssvmCore.calculate_all_scores();
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
