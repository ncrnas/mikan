#include "mk_core_main.hpp"      // MKCoreMain
#include "mke_core.hpp"          // MKECore

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKCoreMain<mkens::MKEOptions, mkens::MKECore>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
