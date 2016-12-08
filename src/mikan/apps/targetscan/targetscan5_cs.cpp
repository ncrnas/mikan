#include <ts5_core.hpp>          // TS5CoreMain

int main(int argc, char const ** argv)
{
    int retVal;

    retVal = ts5cs::TS5CoreMain(argc, argv);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
