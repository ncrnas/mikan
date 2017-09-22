#include "tm1_option.hpp"        // TM1Options

using namespace seqan;

namespace tm1p {

//
// TM1Options methods
//
void TM1Options::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate TargetMiner scores.");
    setVersion(parser, toCString(mProgVer));
    setDate(parser, toCString(mProgDate));

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
                         "\"\\fIOUT_SCORE1 FILE\\fP\" \"\\fIOUT_SCORE2 FILE\\fP\"");
    addDescription(parser,
                   "This program calculates TargetMiner scores for miRNA:mRNA pairs.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "TargetMiner1 Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));
    addOption(parser, ArgParseOption("", "no_gff", "Change output format to tool specific instead of GFF."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBtargetminer1\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_file1\\fP \\fIoutput_file2\\fP",
                "calculate TargetMinder scores of miRNA:mRNA pairs specified in \\fImiRNAs\\fP and \\fImRNA\\fP "
                        "and write the site data to \\fIoutput1\\fP and the scores to \\fIoutput2\\fP");
}


} // namespace tm1p
