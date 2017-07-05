#ifndef TSSVM_OPTION_HPP_
#define TSSVM_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

namespace tssvm{

//
// Tool options
//
class TSSVMOptions
{
public:
    // Define types
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Declare variables
    bool mOutputAlign;
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;
    seqan::CharString mOFileTargetSite;
    seqan::CharString mOFileMRNA;

public:
    // Define methods
    TSSVMOptions() : mOutputAlign(false) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);
    seqan::ArgumentParser::ParseResult validateFiles();
};

} // namespace tssvm

#endif /* TSSVM_OPTION_HPP_ */
