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
    // Declare variables
    bool mOutputAlign;
    int mMinSeedLen;
    int mMaxSeedLen;
    float mMinAlignScore;
    float mMaxEnergy;
    seqan::CharString mAllowGUWobble;
    seqan::CharString mAllowMismatch;
    seqan::CharString mAllowBT;

    // Define methods
    MR3Options() : mOutputAlign(false), mMinSeedLen(6), mMaxSeedLen(8), mMinAlignScore(140.0), mMaxEnergy(1.0) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace mr3as

#endif /* MR3_OPTION_HPP_ */
