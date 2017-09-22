#include "mk_main.hpp"           // MKMain
#include "ts5_core.hpp"          // TS5Core

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKMain<ts5cs::TS5Options, ts5cs::TS5Core>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
