#include "rh2_option.hpp"        // RH2Options

using namespace seqan;

namespace rh2mfe {

//
// RH2Options methods
//
ArgumentParser::ParseResult RH2Options::parseCommandLine(
        int argc,
        char const **argv) {
    // Setup ArgumentParser
    ArgumentParser parser("rnahybrid2_mfe");
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
    mOutputAlign = isSet(parser, "output_align");
    getOptionValue(mSeedDef, parser, "seed_def");
    if (mSeedDef != "6mer" && mSeedDef != "7mer" && mSeedDef != "6mGU+" && mSeedDef != "7mGU+"
        && mSeedDef != "6mGU1" && mSeedDef != "7mGU1") {
        std::cerr << "ERROR: seed_type must be one of the following options: "
                "6mer, 7mer, 6mGU+, 7mGU+, 6mGU1, 7mGU1." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(mOverlapDef, parser, "overlap");
    if (mOverlapDef != "orig" && mOverlapDef != "seed") {
        std::cerr << "ERROR: overlap must be one of the following options: "
                "orig, seed." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(mTargetLen, parser, "target_len");
    getOptionValue(mMaxHits, parser, "max_hits");

    return ArgumentParser::PARSE_OK;
}

void RH2Options::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate RNAhybrid MFE values.");
    setVersion(parser, "1.0");
    setDate(parser, "January 2014");

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
                         "\"\\fIOUT_MFE FILE\\fP\" \"\\fIOUT_TOTAL FILE\\fP\"");
    addDescription(parser,
                   "This program calculates MFE scores and summarizes them for miRNA:mRNA pairs.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "RNAhybrid MFE Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));

    addOption(parser,
              ArgParseOption("s", "seed_def", "Seed site definition [6mer, 7mer, 6mGU+, 7mGU+, 6mGU1, 7mGU1].",
                             ArgParseOption::STRING));
    setDefaultValue(parser, "seed_def", "7mGU+");

    addOption(parser, ArgParseOption("o", "overlap", "Overlapped definition [seed, orig].",
                                     ArgParseOption::STRING));
    setDefaultValue(parser, "overlap", "seed");

    addOption(parser, ArgParseOption("l", "target_len", "Length of input target sequences.",
                                     ArgParseOption::INTEGER));
    setDefaultValue(parser, "target_len", 50);
    setMinValue(parser, "target_len", "25");
//    setMaxValue(parser, "target_len", "100");

    addOption(parser, ArgParseOption("b", "max_hits", "Max number of hits per target (0: Any number) [0].",
                                     ArgParseOption::INTEGER));
    setDefaultValue(parser, "max_hits", 0);
    setMinValue(parser, "max_hits", "0");

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBrnahybrid2_mfe\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_txt1\\fP \\fIoutput_txt2\\fP",
                "calculate MFE scores of \\fImiRNAs\\fP in \\fImRNA\\fP regions "
                        "and write them to \\fIoutput1\\fP and total scores to \\fIoutput2\\fP.");
    addListItem(parser,
                "\\fBrnahybrid2_mfe\\fP \\fB-s 7mer\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_txt1\\fP \\fIoutput_txt2\\fP",
                "calculate MFE scores for targets with at least one 7mer seed site.");
//    addListItem(parser,
//                "\\fBrnahybrid2_mfe\\fP \\fB-no\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
//                "\\fIoutput_txt1\\fP \\fIoutput_txt2\\fP",
//                "select the site with the lowest MEF value when multiple sites are overlapped.");

}

} // namespace rh2mfe
