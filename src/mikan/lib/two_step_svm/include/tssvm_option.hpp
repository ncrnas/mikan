#ifndef TSSVM_OPTION_HPP_
#define TSSVM_OPTION_HPP_

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "mk_option.hpp"        // MKOptions

namespace tssvm {

//
// Tool options
//
class TSSVMOptions : public mikan::MKOptions {
public:
    // Define methods
    TSSVMOptions() {
        mProgName = "mk-twostep-svm";
    }

private:
    void setProgramDescription(seqan::ArgumentParser &pParser);
};

} // namespace tssvm

#endif /* TSSVM_OPTION_HPP_ */
