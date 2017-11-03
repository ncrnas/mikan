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
    mGff = isSet(parser, "gff");
    mShowConfig = isSet(parser, "show_config");

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
                 "[\\fIOPTIONS\\fP] \\fImirna_fasta_file\\fP \\fImrna_fasta_file\\fP "
                         "\\fIsite_output_file\\fP \\fIrna_output_file\\fP");
    addDescription(parser,
                   "This program calculates mikan ensemble scores for candidates of miRNA targets.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "Options");
    addOption(parser, ArgParseOption("", "gff", "Change output format to GFF."));

    addOption(parser, ArgParseOption("c", "config", "Specify a configuration file.",
                                     ArgParseArgument::INPUTFILE));
    setDefaultValue(parser, "config", "");

    addOption(parser, ArgParseOption("", "show_config", "Show the content of the specified configuration file."));
    hideOption(parser, "show_config");

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmikan\\fP \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "read two input files and write mikan ensemble scores in two output files.");
    addListItem(parser,
                "\\fBfBmikan\\fP --gff \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "change output format to GFF.");
    addListItem(parser,
                "\\fBfBmikan\\fP -c \\fImikan.conf\\fP \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "specify a configuration file.");

}


} // namespace mkens

