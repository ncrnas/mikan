#ifndef PITA_OPTION_HPP_
#define PITA_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <mk_option.hpp>        // MKOptions

namespace ptddg {

//
// Tool options
//
class PITAOptions : public mikan::MKOptions {
public:
    // Declare variables
    bool mOutputAlign;
    int mMinSeedLen;
    int mMaxSeedLen;
    int mFlankUp;
    int mFlankDown;
    seqan::CharString mAllowGUWobble;
    seqan::CharString mAllowMismatch;

    // Define methods
    PITAOptions() : mOutputAlign(false), mMinSeedLen(6), mMaxSeedLen(8), mFlankUp(0), mFlankDown(0) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace ptddg

#endif /* PITA_OPTION_HPP_ */
