#ifndef MR3_OPTION_HPP_
#define MR3_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

namespace mr3as{

//
// Tool options
//
class MR3Options
{
public:
    // Define types
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Declare variables
    bool mOutputAlign;
    int mMinSeedLen;
    int mMaxSeedLen;
    float mMinAlignScore;
    float mMaxEnergy;
    seqan::CharString mAllowGUWobble;
    seqan::CharString mAllowMismatch;
    seqan::CharString mAllowBT;
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;
    seqan::CharString mOFileSite;
    seqan::CharString mOFileTotal;

public:
    // Define methods
    MR3Options() : mOutputAlign (false), mMinSeedLen(6), mMaxSeedLen(8), mMinAlignScore(140.0), mMaxEnergy(1.0) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);
    seqan::ArgumentParser::ParseResult validateFiles();
};

} // namespace mr3as

#endif /* MR3_OPTION_HPP_ */
