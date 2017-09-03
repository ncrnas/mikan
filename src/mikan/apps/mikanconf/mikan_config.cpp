#include "mke_config.hpp"          // MKEConfig

int main(int, char const **) {
    mkens::MKEConfig conf;

    conf.print_config();

    return 0;
}
