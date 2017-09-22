#include "mk_main.hpp"           // MKMain
#include "pita_core.hpp"         // PITACore

int main(int argc, char const **argv) {
    int retVal;

    retVal = mikan::MKMain<ptddg::PITAOptions, ptddg::PITACore>(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    return 0;
}
