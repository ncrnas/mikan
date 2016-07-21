#include <mikan/lib/rnahybrid2_mfe/include/rh2_core.hpp>          // RH2CoreInput, RH2Core
#include <mikan/lib/rnahybrid2_mfe/include/rh2_inst_template.hpp> // TRNATYPE
#include <mikan/lib/rnahybrid2_mfe/include/rh2_option.hpp>        // RH2Options
#include <seqan/arg_parse.h>

int main(int argc, char const ** argv)
{
    int retVal;

    // Parse the command line.
    rh2mfe::RH2Options options;
    seqan::ArgumentParser::ParseResult parseRes = options.parseCommandLine(argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
    {
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read input files
    rh2mfe::RH2CoreInput<rh2mfe::TRNATYPE> coreInput;
    coreInput.init_from_args(options);
    retVal = coreInput.load_seq_from_file();
    if (retVal != 0)
    {
        return 1;
    }

    // Create index
    rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TRNASet const& mMRNASeqs = coreInput.get_mrna_seqs();
    rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TIndexQGram index(mMRNASeqs);
    rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TFinder finder(index);

    // Calculate scores for all miRNAs
    rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TCharSet const& mMiRNAIds = coreInput.get_mirna_ids();
    rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TRNASet const& mMiRNASeqs = coreInput.get_mirna_seqs();
    rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TCharSet const& mMRNAIds = coreInput.get_mrna_ids();
    int mRNAMaxLen = options.mTargetLen;
    int miRNAMaxLen = options.mQueryLen;
    std::string seedDef(toCString(options.mSeedDef));
    rh2mfe::RH2Core<rh2mfe::TRNATYPE> rh2Core(mMiRNAIds, mMiRNASeqs, mMRNAIds, mMRNASeqs, index, finder,
            mRNAMaxLen, miRNAMaxLen, seedDef);
    rh2Core.init_from_args(options);
    rh2Core.open_output_file();
    retVal = rh2Core.calculate_all_scores();
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
