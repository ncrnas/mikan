#include "mk_main.hpp"           // MKMain
#include "rh2_core.hpp"          // RH2Core

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKMain<rh2mfe::RH2Options, rh2mfe::RH2Core>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
