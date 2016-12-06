#include <pita_core.hpp>          // PITACoreMain

int main(int argc, char const ** argv)
{
    int retVal;

    retVal = ptddg::PITACoreMain(argc, argv);
    if (retVal != 0)
    {
        return 1;
    }

    return 0;
}
