#include "mk_core_main.hpp"      // MKCoreMain
#include "rh2_core.hpp"          // RH2Core

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKCoreMain<rh2mfe::RH2Options, rh2mfe::RH2Core>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
