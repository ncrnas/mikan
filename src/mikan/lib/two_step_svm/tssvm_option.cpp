#include "tssvm_option.hpp"         // TSSVMOptions

using namespace seqan;

namespace tssvm {

//
// TSSVMOptions methods
//
void TSSVMOptions::setProgramDescription(seqan::ArgumentParser &parser) {
    // Set short description, version, and date
    setShortDescription(parser, "Calculate Two-step SVM scores.");
    setVersion(parser, toCString(mProgVer));
    setDate(parser, toCString(mProgDate));

    // Define usage line and long description
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fImirna_fasta_file\\fP \\fImrna_fasta_file\\fP "
                         "\\fIsite_output_file\\fP \\fIrna_output_file\\fP");
    addDescription(parser,
                   "This program calculates Two-step SVM scores for candidates of miRNA targets.");

    // Define Arguments
    addIOArgs(parser);

    // Define Options
    addSection(parser, "Options");
    addOption(parser, ArgParseOption("a", "output_align", "Output alignments to standard output."));
    addOption(parser, ArgParseOption("", "gff", "Change output format to GFF."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBmk-twostep-svm\\fP \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "read two input files and write Two-step SVM scores in two output files.");
    addListItem(parser,
                "\\fBmk-twostep-svm\\fP -a \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "output alignments to standard output in addition.");
    addListItem(parser,
                "\\fBmk-twostep-svm\\fP --gff \\fImirna.fasta\\fP \\fImrna.fasta\\fP "
                        "\\fIsite_output.txt\\fP \\fIrna_output.txt\\fP",
                "change output format to GFF.");
}

} // namespace tssvm
