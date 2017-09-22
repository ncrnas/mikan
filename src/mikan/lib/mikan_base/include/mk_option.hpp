#ifndef MK_OPTION_HPP_
#define MK_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_const.hpp"
#include "mk_typedef.hpp"

namespace mikan {

//
// Tool options
//
class MKOptions {
public:
    // Define type
    typedef seqan::ArgumentParser::ParseResult TParseResult;

    // Define variables
    mikan::TCharStr mMiRNAFasta;
    mikan::TCharStr mMRNAFasta;
    mikan::TCharStr mOFileSite;
    mikan::TCharStr mOFileTotal;

    bool mOutputAlign;
    bool mNoGff;

    // MR3 & PITA
    int mMinSeedLen;
    int mMaxSeedLen;
    mikan::TCharStr mAllowGUWobble;
    mikan::TCharStr mAllowMismatch;

    // MR3
    float mMinAlignScore;
    float mMaxEnergy;
    mikan::TCharStr mAllowBT;

    // PITA
    int mFlankUp;
    int mFlankDown;

    // RH2
    std::string mSeedDef;
    mikan::TCharStr mOverlapDef;
    int mTargetLen;
    int mQueryLen;
    int mMaxHits;

    // MKRMAWithSites
    mikan::TCharStr mSortValType;

    // Mikan ensemble
    mikan::TCharSet mToolPrefix;
    mikan::TCharStr mConfigFile;
    bool mShowConfig;

    // Define method
    MKOptions() {
        mProgName = "mikan";
        mProgVer = "1.0";
        mProgDate = "July 2017";

        mOutputAlign = false;
        mNoGff = false;

        mMaxHits = 0;

        mSortValType = "position";

        mConfigFile = "";
        mShowConfig = false;
    }

    // Method prototype
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

protected:
    // Define variables
    mikan::TCharStr mProgName;
    mikan::TCharStr mProgVer;
    mikan::TCharStr mProgDate;

    // Method prototypes
    seqan::ArgumentParser::ParseResult validateFiles(seqan::ArgumentParser &parser);

    static void addIOArgs(seqan::ArgumentParser &parser);

private:
    // Method prototype
    virtual void setProgramDescription(seqan::ArgumentParser &pParser);

};

} // namespace mikan

#endif /* MK_OPTION_HPP_ */
