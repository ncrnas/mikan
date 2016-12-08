#include <tm1_core.hpp>          // TM1CoreMain

int main(int argc, char const ** argv)
{
    int retVal;

    retVal = tm1p::TM1CoreMain(argc, argv);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
