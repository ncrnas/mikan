#ifndef TM1_OPTION_HPP_
#define TM1_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

namespace tm1p {

//
// Tool options
//
class TM1CSOptions {
public:
    // Define types
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Declare variables
    bool mOutputAlign;
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileScore;

public:
    // Define methods
    TM1CSOptions() : mOutputAlign(false) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);

    seqan::ArgumentParser::ParseResult validateFiles();
};

} // namespace tm1p

#endif /* TM1_OPTION_HPP_ */
