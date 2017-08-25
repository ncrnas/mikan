#ifndef MK_OPTION_HPP_
#define MK_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_const.hpp"

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

    bool mOutputAlign;

    // MR3 & PITA
    int mMinSeedLen;
    int mMaxSeedLen;
    seqan::CharString mAllowGUWobble;
    seqan::CharString mAllowMismatch;

    // MR3
    float mMinAlignScore;
    float mMaxEnergy;
    seqan::CharString mAllowBT;

    // PITA
    int mFlankUp;
    int mFlankDown;

    // RH2
    std::string mSeedDef;
    seqan::CharString mOverlapDef;
    int mTargetLen;
    int mQueryLen;
    int mMaxHits;

    // MKRMAWithSites
    seqan::CharString mSortValType;

    // Mikan ensemble
    seqan::StringSet<seqan::CharString> mToolPrefix;

    // Define method
    MKOptions() {
        mProgName = "mikan";
        mProgVer = "1.0";
        mProgDate = "July 2017";

        mOutputAlign = false;

        mMaxHits = 0;

        mSortValType = "position";
    }

    // Method prototype
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

protected:
    // Define variables
    seqan::CharString mProgName;
    seqan::CharString mProgVer;
    seqan::CharString mProgDate;

    // Method prototypes
    seqan::ArgumentParser::ParseResult validateFiles(seqan::ArgumentParser &parser);

    static void addIOArgs(seqan::ArgumentParser &parser);

private:
    // Method prototype
    virtual void setProgramDescription(seqan::ArgumentParser &pParser);

};

} // namespace mikan

#endif /* MK_OPTION_HPP_ */
