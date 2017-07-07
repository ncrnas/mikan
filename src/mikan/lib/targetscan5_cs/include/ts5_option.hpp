#ifndef TS5_OPTION_HPP_
#define TS5_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions

namespace ts5cs {

//
// Tool options
//
class TS5CSOptions : public mikan::MKOptions {
public:
    // Declare variables
    bool mOutputAlign;

    // Define methods
    TS5CSOptions() : mOutputAlign(false) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace ts5cs

#endif /* TS5_OPTION_HPP_ */
