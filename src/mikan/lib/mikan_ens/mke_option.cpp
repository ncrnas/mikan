#include "mke_option.hpp"        // MKEOptions

using namespace seqan;

namespace mkens {

//
// MKEOptions methods
//
ArgumentParser::ParseResult MKEOptions::parseCommandLine(
        int argc,
        char const **argv) {

    // Setup ArgumentParser
    ArgumentParser parser(toCString(mProgName));
    setProgramDescription(parser);

    // Parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK) {
        return res;
    }

    // Validate files
    res = validateFiles(parser);
    if (res != ArgumentParser::PARSE_OK) {
        return res;
    }

    // Extract options
    mShowConfig = isSet(parser, "show-config");

    getOptionValue(mConfigFile, parser, "config");
    if (mConfigFile != "") {
        char const *conf = toCString(mConfigFile);
        std::fstream conf_file(conf, std::ios::in);
        if (!conf_file.good()) {
            std::cerr << "ERROR: Could not open the specified configuration file: " << conf << std::endl;
            return ArgumentParser::PARSE_ERROR;
        }

        std::string conf_f = conf;
        mConf.parse_config(conf_f);
        if (mShowConfig) {
            mConf.print_config();
        }
    }

    return ArgumentParser::PARSE_OK;
}

void MKEOptions::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate mikan ensemble scores.");
    setVersion(parser, toCString(mProgVer));
    setDate(parser, toCString(mProgDate));

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
                         "\"\\fIOUT_SCORE1 FILE\\fP\" \"\\fIOUT_SCORE2 FILE\\fP\"");
    addDescription(parser,
                   "This program calculates mikan ensemble scores and summarizes them for miRNA:mRNA pairs.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "Mikan Ensemble Score Options");
    addOption(parser, ArgParseOption("c", "config", "Specify a configuration file.",
                                     ArgParseArgument::INPUTFILE));
    setDefaultValue(parser, "config", "");

    addOption(parser, ArgParseOption("", "show-config", "Show the content of the specified configuration file."));
    hideOption(parser, "show-config");

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmikan\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_file1\\fP \\fIoutput_file2\\fP",
                "calculate context scores of \\fImiRNAs\\fP in \\fImRNA\\fP regions "
                        "and write the scores to \\fIoutput1\\fP and the total scores to \\fIoutput2\\fP");
}


} // namespace mkens

