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
    // Define type
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Define variables
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileTotal;

    // Define method
    explicit MKOptions() {}

    // Method prototype
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

protected:
    // Method prototypes
    seqan::ArgumentParser::ParseResult validateFiles(seqan::ArgumentParser &parser);

    static void addIOArgs(seqan::ArgumentParser &parser);

private:
    // Method prototype
    static void setProgramDescription(seqan::ArgumentParser &pParser);

};

} // namespace mikan

#endif /* MK_OPTION_HPP_ */
