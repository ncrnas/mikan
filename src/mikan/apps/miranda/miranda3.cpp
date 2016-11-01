#include <mr3_core.hpp>          // MR3CoreMain

int main(int argc, char const ** argv)
{
    int retVal;

    retVal = mr3as::MR3CoreMain(argc, argv);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
