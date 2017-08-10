#ifndef MKE_OPTION_HPP_
#define MKE_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions

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

private:
    virtual void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace mkens

#endif /* MKE_OPTION_HPP_ */
