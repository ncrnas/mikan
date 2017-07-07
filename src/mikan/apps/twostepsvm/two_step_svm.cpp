#include "tssvm_core.hpp"           // TSSVMCoreMain

int main(int argc, char const **argv) {
    int retVal;

    retVal = tssvm::TSSVMCoreMain(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
