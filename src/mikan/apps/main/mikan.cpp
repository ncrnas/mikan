#include "mr3_core.hpp"          // MR3CoreMain
#include "pita_core.hpp"         // PITACoreMain
#include "rh2_core.hpp"          // RH2CoreMain
#include "tm1_core.hpp"          // TM1CoreMain
#include "ts5_core.hpp"          // TS5CoreMain
#include "tssvm_core.hpp"        // TSSVMCoreMain

int main(int argc, char const **argv) {
    int retVal;

    retVal = mr3as::MR3CoreMain(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    retVal = ptddg::PITACoreMain(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    retVal = rh2mfe::RH2CoreMain(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    retVal = tm1p::TM1CoreMain(argc, argv);
    if (retVal != 0) {
        return 1;
    }

    retVal = ts5cs::TS5CoreMain(argc, argv);
    if (retVal != 0) {
        return 1;
    }


    return 0;
}
