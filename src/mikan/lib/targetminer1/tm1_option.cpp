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
                 "[\\fIOPTIONS\\fP] \\fImirna_fasta_file\\fP \\fImrna_fasta_file\\fP "
                         "\\fIsite_output_file\\fP \\fIrna_output_file\\fP");
    addDescription(parser,
                   "This program calculates TargetMiner scores for candidates of miRNA targets.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));
    addOption(parser, ArgParseOption("", "gff", "Change output format to GFF."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmk-targetminer\\fP \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "read two input files and write TargetMiner scores in two output files.");
    addListItem(parser,
                "\\fBmk-targetminer\\fP -a \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "output alignments to standard output in addition.");
    addListItem(parser,
                "\\fBmk-targetminer\\fP --gff \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "change output format to GFF.");
}


} // namespace tm1p
