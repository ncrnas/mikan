#ifndef MK_OPTION_HPP_
#define MK_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

namespace mikan {

//
// Tool options
//
class MKOptions {
public:
    // Define types
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Declare variables
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileTotal;

public:
    // Define methods
    MKOptions() {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);

protected:
    seqan::ArgumentParser::ParseResult validateFiles(seqan::ArgumentParser &parser);

    static void addIOArgs(seqan::ArgumentParser &parser);
};

} // namespace mikan

#endif /* MK_OPTION_HPP_ */
