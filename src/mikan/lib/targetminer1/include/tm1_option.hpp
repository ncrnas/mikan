#ifndef TM1_OPTION_HPP_
#define TM1_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions

namespace tm1p {

//
// Tool options
//
class TM1CSOptions : public mikan::MKOptions {
public:
    // Define methods
    TM1CSOptions() {
        mProgName = "mk-targetminer";
    }

private:
    void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace tm1p

#endif /* TM1_OPTION_HPP_ */
