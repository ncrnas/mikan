#include <mikan/lib/targetminer1/include/tm1_option.hpp>        // TM1CSOptions

using namespace seqan;

namespace tm1p{

//
// TM1CSOptions methods
//
ArgumentParser::ParseResult TM1CSOptions::parseCommandLine(
        int argc,
        char const **argv)
{
    // Setup ArgumentParser
    ArgumentParser parser("targetminer1");
    setProgramDescription(parser);

    // Parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        return res;
    }

    // Extract arguments
    getArgumentValue(mMiRNAFasta, parser, 0);
    getArgumentValue(mMRNAFasta, parser, 1);
    getArgumentValue(mOFileSite, parser, 2);
    getArgumentValue(mOFileScore, parser, 3);

    // Validate files
    res = validateFiles();
    if (res != ArgumentParser::PARSE_OK)
    {
        return res;
    }

    // Extract options
    mOutputAlign = isSet(parser, "output_align");

    return ArgumentParser::PARSE_OK;
}

void TM1CSOptions::setProgramDescription(seqan::ArgumentParser &parser)
{
    // Set short description, version, and date
    setShortDescription(parser, "Calculate TargetMiner scores.");
    setVersion(parser, "1.0");
    setDate(parser, "January 2014");

    // Define usage line and long description
    addUsageLine(parser,
            "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
            "\"\\fIOUT_SCORE1 FILE\\fP\" \"\\fIOUT_SCORE2 FILE\\fP\"");
    addDescription(parser,
            "This program calculates TargetMiner scores for miRNA:mRNA pairs.");

    // Define Arguments
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE));

    // Define Options
    addSection(parser, "TargetMiner1 Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
            "\\fBtargetminer1\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
            "\\fIoutput_file1\\fP \\fIoutput_file2\\fP",
            "calculate TargetMinder scores of miRNA:mRNA pairs specified in \\fImiRNAs\\fP and \\fImRNA\\fP "
            "and write the site data to \\fIoutput1\\fP and the scores to \\fIoutput2\\fP");
}

ArgumentParser::ParseResult TM1CSOptions::validateFiles()
{
    char const *input_1 = toCString(mMiRNAFasta);
    char const *input_2 = toCString(mMRNAFasta);
    char const *output_1 = toCString(mOFileSite);
    char const *output_2 = toCString(mOFileScore);

    std::fstream mirna_fa(input_1, std::ios::in);
    if (!mirna_fa.good())
    {
        std::cerr << "ERROR: Could not open input file: " << input_1 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    std::fstream mrna_fa(input_2, std::ios::in);
    if (!mrna_fa.good())
    {
        std::cerr << "ERROR: Could not open input file: " << input_2 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    std::fstream ofile1(output_1, std::ios::out);
    if (!ofile1.good())
    {
        std::cerr << "ERROR: Could not open output file: " << output_1 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    std::fstream ofile2(output_2, std::ios::out);
    if (!ofile2.good())
    {
        std::cerr << "ERROR: Could not open output file: " << output_2 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    return ArgumentParser::PARSE_OK;
}

} // namespace tm1p
