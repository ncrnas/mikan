#ifndef MIKAN_TEST_PITA_HPP_
#define MIKAN_TEST_PITA_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mk_main.hpp"
#include "pita_core.hpp"
#include "mk_input.hpp"

typedef TestSeed<ptddg::PITAOptions, ptddg::PITASeedSeqs> TestSeedPITA;
typedef TestSite<ptddg::PITAOptions, ptddg::PITASeedSeqs, ptddg::PITASeedSites> TestSitePITA;

#endif //MIKAN_TEST_PITA_HPP_
