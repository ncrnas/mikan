#ifndef TS5_OPTION_HPP_
#define TS5_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

namespace ts5cs{

//
// Tool options
//
class TS5CSOptions
{
public:
    // Define types
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Declare variables
    bool mOutputAlign;
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;
    seqan::CharString mOFileContext;
    seqan::CharString mOFileTotal;

public:
    // Define methods
    TS5CSOptions() : mOutputAlign(false) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);
    seqan::ArgumentParser::ParseResult validateFiles();
};

} // namespace ts5cs

#endif /* TS5_OPTION_HPP_ */
