#include <mk_option.hpp>         // MKOptions

using namespace seqan;

namespace mikan {

//
// MKOptions methods
//
ArgumentParser::ParseResult MKOptions::parseCommandLine(
        int argc,
        char const **argv) {
    // Setup ArgumentParser
    ArgumentParser parser("mikan");
    setProgramDescription(parser);

    // Parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK) {
        return res;
    }

    return ArgumentParser::PARSE_OK;
}

void MKOptions::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate mikan scores.");
    setVersion(parser, "1.0");
    setDate(parser, "July 2017");

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
                         "\"\\fIOUT_DDG FILE\\fP\" \"\\fIOUT_TOTAL FILE\\fP\"");
    addDescription(parser,
                   "This program calculates mikan scores and summarizes them for miRNA:mRNA pairs.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
//    addSection(parser, "mikan Options");

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmikan\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_txt1\\fP \\fIoutput_txt2\\fP",
                "calculate mikan scores of \\fImiRNAs\\fP in \\fImRNA\\fP regions "
                        "and write them to \\fIoutput1\\fP and total scores to \\fIoutput2\\fP.");

}

ArgumentParser::ParseResult MKOptions::validateFiles(seqan::ArgumentParser &parser) {
    // Extract arguments
    getArgumentValue(mMiRNAFasta, parser, 0);
    getArgumentValue(mMRNAFasta, parser, 1);
    getArgumentValue(mOFileSite, parser, 2);
    getArgumentValue(mOFileTotal, parser, 3);

    char const *input_1 = toCString(mMiRNAFasta);
    char const *input_2 = toCString(mMRNAFasta);
    char const *output_1 = toCString(mOFileSite);
    char const *output_2 = toCString(mOFileTotal);

    std::fstream mirna_fa(input_1, std::ios::in);
    if (!mirna_fa.good()) {
        std::cerr << "ERROR: Could not open input file: " << input_1 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    std::fstream mrna_fa(input_2, std::ios::in);
    if (!mrna_fa.good()) {
        std::cerr << "ERROR: Could not open input file: " << input_2 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    std::fstream ofile1(output_1, std::ios::out);
    if (!ofile1.good()) {
        std::cerr << "ERROR: Could not open output file: " << output_1 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    std::fstream ofile2(output_2, std::ios::out);
    if (!ofile2.good()) {
        std::cerr << "ERROR: Could not open output file: " << output_2 << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    return ArgumentParser::PARSE_OK;
}

void MKOptions::addIOArgs(seqan::ArgumentParser &parser) {
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE));
}

} // namespace mikan
