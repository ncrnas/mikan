#include "mk_core_main.hpp"      // MKCoreMain
#include "tssvm_core.hpp"        // TSSVMCore

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKCoreMain<tssvm::TSSVMOptions, tssvm::TSSVMCore>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
