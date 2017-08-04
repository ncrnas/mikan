#ifndef MR3_OPTION_HPP_
#define MR3_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions

namespace mr3as {

//
// Tool options
//
class MR3Options : public mikan::MKOptions {
public:
    // Define method
    MR3Options() {
        mProgName = "mk-mirnada";

        mOutputAlign = 6;
        mMaxSeedLen = 8;
        mMinAlignScore = 140.0;
        mMaxEnergy = 1.0;
    }

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    virtual void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace mr3as

#endif /* MR3_OPTION_HPP_ */
