#include "mk_core_main.hpp"      // MKCoreMain
#include "pita_core.hpp"         // PITACore

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKCoreMain<ptddg::PITAOptions, ptddg::PITACore >(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
