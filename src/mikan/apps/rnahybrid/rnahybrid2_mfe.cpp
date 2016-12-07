#include <rh2_core.hpp>          // RH2CoreMain

int main(int argc, char const ** argv)
{
    int retVal;

    retVal = rh2mfe::RH2CoreMain(argc, argv);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
