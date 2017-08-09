#include "mk_core_main.hpp"      // MKCoreMain
#include "mr3_core.hpp"          // MR3Core

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKCoreMain<mr3as::MR3Options, mr3as::MR3Core >(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
