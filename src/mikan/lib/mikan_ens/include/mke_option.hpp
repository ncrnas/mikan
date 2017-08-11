#ifndef MKE_OPTION_HPP_
#define MKE_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions
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
    }

    mr3as::MR3Options mMR3Opts;
    ptddg::PITAOptions mPITAOpts;
    rh2mfe::RH2Options mRH2Opts;
    tm1p::TM1Options mTM1Opts;
    ts5cs::TS5Options mTS5Opts;
    tssvm::TSSVMOptions mTSSVMOpts;

private:
    virtual void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace mkens

#endif /* MKE_OPTION_HPP_ */
