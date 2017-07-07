#include "mr3_option.hpp"         // MR3Options

using namespace seqan;

namespace mr3as {

//
// MR3Options methods
//
ArgumentParser::ParseResult MR3Options::parseCommandLine(
        int argc,
        char const **argv) {
    // Setup ArgumentParser
    ArgumentParser parser("mirnada3");
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

    getOptionValue(mMinSeedLen, parser, "seed_min");
    if (mMinSeedLen != 6 && mMinSeedLen != 7 && mMinSeedLen != 8) {
        std::cerr << "ERROR: seed_min must be one of the following options: 6, 7, 8" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(mMaxSeedLen, parser, "seed_max");
    if (mMaxSeedLen != 6 && mMaxSeedLen != 7 && mMaxSeedLen != 8) {
        std::cerr << "ERROR: seed_max must be one of the following options: 6, 7, 8" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    if (mMinSeedLen > mMaxSeedLen) {
        std::cerr << "ERROR: seed_min must be equal to or smaller than seed_max." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(mAllowGUWobble, parser, "allow_gu");
    if (mAllowGUWobble != "0" && mAllowGUWobble != "1" && mAllowGUWobble != "+") {
        std::cerr << "ERROR: allow_gu must be one of the following options: 0, 1, +" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(mAllowMismatch, parser, "allow_mm");
    if (mAllowMismatch != "0:0" && mAllowMismatch != "0:1" && mAllowMismatch != "1:1") {
        std::cerr << "ERROR: allow_mm must be one of the following options: 0:0, 0:1, 1:1" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(mAllowBT, parser, "allow_bt");
    if (mAllowBT != "0" && mAllowBT != "1") {
        std::cerr << "ERROR: allow_bt must be one of the following options: 0, 1" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    return ArgumentParser::PARSE_OK;
}

void MR3Options::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate miRanda scores.");
    setVersion(parser, "1.0");
    setDate(parser, "January 2014");

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
                         "\"\\fIOUT_DDG FILE\\fP\" \"\\fIOUT_TOTAL FILE\\fP\"");
    addDescription(parser,
                   "This program calculates miRanda scores and summarizes them for miRNA:mRNA pairs.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "miRanda3 Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));

    addOption(parser,
              ArgParseOption("n", "seed_min", "Seed site minimum length [6, 7, 8].", ArgParseOption::INTEGER));
    setDefaultValue(parser, "seed_min", 6);

    addOption(parser,
              ArgParseOption("x", "seed_max", "Seed site maximum length [6, 7, 8].", ArgParseOption::INTEGER));
    setDefaultValue(parser, "seed_max", 8);

    addOption(parser, ArgParseOption("l", "min_align", "Minimum alignment score.", ArgParseOption::INTEGER));
    setDefaultValue(parser, "min_align", 140);

    addOption(parser, ArgParseOption("e", "max_en", "Maximum energy score.", ArgParseOption::DOUBLE));
    setDefaultValue(parser, "max_en", 1);

    addOption(parser, ArgParseOption("g", "allow_gu", "Allow single(1) or multiple(+) GU wobbles for 7 or 8mers "
            "[0, 1, +].", ArgParseOption::STRING));
    setDefaultValue(parser, "allow_gu", "+");

    addOption(parser, ArgParseOption("m", "allow_mm", "Allow one mismatch for 7 or 8mers [0:0, 0:1, 1:1].",
                                     ArgParseOption::STRING));
    setDefaultValue(parser, "allow_mm", "1:1");

    addOption(parser, ArgParseOption("t", "allow_bt", "Allow one bulge in position 2-7 of target site [0, 1].",
                                     ArgParseOption::STRING));
    setDefaultValue(parser, "allow_bt", "1");

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmirnada3\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_txt1\\fP \\fIoutput_txt2\\fP",
                "calculate miRanda alignment and energy scores of \\fImiRNAs\\fP in \\fImRNA\\fP regions "
                        "and write them to \\fIoutput1\\fP and total scores to \\fIoutput2\\fP.");

}

} // namespace mr3as
