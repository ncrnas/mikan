#include <tssvm_option.hpp>         // TSSVMOptions

using namespace seqan;

namespace tssvm {

//
// TSSVMOptions methods
//
ArgumentParser::ParseResult TSSVMOptions::parseCommandLine(
        int argc,
        char const **argv) {
    // Setup ArgumentParser
    ArgumentParser parser("two_step_svm");
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

    return ArgumentParser::PARSE_OK;
}

void TSSVMOptions::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate Two-step SVM scores.");
    setVersion(parser, "1.0");
    setDate(parser, "January 2014");

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
                         "\"\\fIOUT_TARGETSITE FILE\\fP\" \"\\fIOUT_MRNA\\fP\"");
    addDescription(parser,
                   "This program calculates Two-step SVM scores for both target-site and mRNA levels.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "Two-step SVM Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBtwo_step_svm\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_file1\\fP \\fIoutput_file2\\fP",
                "calculate Two-step SVM scores of \\fImiRNAs\\fP in \\fImRNA\\fP regions "
                        "by using models in \\fImodel_path\\fP folder "
                        "and write the target-site level scores to \\fIoutput1\\fP "
                        "and mRNA level scores to \\fIoutput2\\fP");
}

} // namespace tssvm
