#include "ts5_option.hpp"        // TS5Options

using namespace seqan;

namespace ts5cs {

//
// TS5Options methods
//
void TS5Options::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate TargetScan 5 context scores.");
    setVersion(parser, toCString(mProgVer));
    setDate(parser, toCString(mProgDate));

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" "
                         "\"\\fIOUT_SCORE1 FILE\\fP\" \"\\fIOUT_SCORE2 FILE\\fP\"");
    addDescription(parser,
                   "This program calculates TargetScan context scores and summarizes them for miRNA:mRNA pairs.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "TargetScan5 Context Score Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));
    addOption(parser, ArgParseOption("", "no_gff", "Change output format to tool specific instead of GFF."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBtargetscan5_cs\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_file1\\fP \\fIoutput_file2\\fP",
                "calculate context scores of \\fImiRNAs\\fP in \\fImRNA\\fP regions "
                        "and write the scores to \\fIoutput1\\fP and the total scores to \\fIoutput2\\fP");
}


} // namespace ts5cs
