#include "mk_core_main.hpp"      // MKCoreMain
#include "tm1_core.hpp"          // TM1Core

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKCoreMain<tm1p::TM1Options, tm1p::TM1Core >(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
