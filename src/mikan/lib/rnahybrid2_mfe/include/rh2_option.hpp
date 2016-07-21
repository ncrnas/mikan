#ifndef RH2_OPTION_HPP_
#define RH2_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

namespace rh2mfe{

//
// Tool options
//
class RH2Options
{
public:
    // Define types
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Declare variables
    bool mOutputAlign;
    seqan::CharString mSeedDef;
    seqan::CharString mOverlapDef;
    int mTargetLen;
    int mQueryLen;
    seqan::CharString mMiRNAFasta;
    seqan::CharString mMRNAFasta;
    seqan::CharString mOFileMFE;
    seqan::CharString mOFileTotal;
    int mMaxHits;

public:
    // Define methods
    RH2Options() : mOutputAlign (false), mSeedDef("7mGU+"), mOverlapDef("orig"), mTargetLen(50), mQueryLen(30),
    mMaxHits(0) {}

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    static void setProgramDescription(seqan::ArgumentParser &pParser);
    seqan::ArgumentParser::ParseResult validateFiles();
};

} // namespace rh2mfe

#endif /* RH2_OPTION_HPP_ */
