#ifndef MKE_OPTION_HPP_
#define MKE_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions
#include "mke_const.hpp"        // TOOL_NUM
#include "mke_config.hpp"       // MKEConfig
#include "mr3_option.hpp"       // MR3Options
#include "pita_option.hpp"      // PITAOptions
#include "rh2_option.hpp"       // RH2Options
#include "tm1_option.hpp"       // TM1Options
#include "ts5_option.hpp"       // TS5Options
#include "tssvm_option.hpp"     // TSSVMOptions

namespace mkens {

//
// Tool options
//
class MKEOptions : public mikan::MKOptions {
public:
    // Define method
    MKEOptions() {
        mProgName = "mikan";

        resize(mToolPrefix, mkens::TOOL_NUM);
        mToolPrefix[0] = "mr";
        mToolPrefix[1] = "pt";
        mToolPrefix[2] = "rh";
        mToolPrefix[3] = "tm";
        mToolPrefix[4] = "ts";
        mToolPrefix[5] = "sv";

    }

    mr3as::MR3Options mMR3Opts;
    ptddg::PITAOptions mPITAOpts;
    rh2mfe::RH2Options mRH2Opts;
    tm1p::TM1Options mTM1Opts;
    ts5cs::TS5Options mTS5Opts;
    tssvm::TSSVMOptions mTSSVMOpts;

    MKEConfig mConf;

    // Method prototype
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const **argv);

private:
    virtual void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace mkens

#endif /* MKE_OPTION_HPP_ */
