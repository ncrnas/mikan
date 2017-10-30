#include "mk_main.hpp"           // MKMain
#include "mke_core.hpp"          // MKECore

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKMain<mkens::MKEOptions, mkens::MKECore>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
