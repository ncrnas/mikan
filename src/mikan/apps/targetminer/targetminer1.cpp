#include "mk_main.hpp"           // MKMain
#include "tm1_core.hpp"          // TM1Core

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
