#ifndef TS5_OPTION_HPP_
#define TS5_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions

namespace ts5cs {

//
// Tool options
//
class TS5CSOptions : public mikan::MKOptions {
public:
    // Define method
    TS5CSOptions() {
        mProgName = "mk-targetscan";
    }

private:
    virtual void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace ts5cs

#endif /* TS5_OPTION_HPP_ */
