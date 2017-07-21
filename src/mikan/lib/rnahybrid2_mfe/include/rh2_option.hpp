#ifndef RH2_OPTION_HPP_
#define RH2_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions

namespace rh2mfe {

//
// Tool options
//
class RH2Options : public mikan::MKOptions {
public:
    // Define method
    RH2Options() {
        mProgName = "mk-rnahybrid";

        mSeedDef = "7mGU+";
        mOverlapDef = "seed";
        mTargetLen = 50;
        mQueryLen = 30;
        mMaxHits = 0;
    }

    // Method prototypes
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace rh2mfe

#endif /* RH2_OPTION_HPP_ */
