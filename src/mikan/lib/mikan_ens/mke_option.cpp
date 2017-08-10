#include "mke_option.hpp"        // MKEOptions

using namespace seqan;

namespace mkens {

//
// MKEOptions methods
//
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
    addSection(parser, "mikan ensemble Score Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmikan\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP "
                        "\\fIoutput_file1\\fP \\fIoutput_file2\\fP",
                "calculate context scores of \\fImiRNAs\\fP in \\fImRNA\\fP regions "
                        "and write the scores to \\fIoutput1\\fP and the total scores to \\fIoutput2\\fP");
}


} // namespace mkens

