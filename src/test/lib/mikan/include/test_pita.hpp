#ifndef MIKAN_TEST_PITA_HPP_
#define MIKAN_TEST_PITA_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mk_core_main.hpp"
#include "pita_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput> TestIOPITA;
typedef TestSeed<ptddg::PITASeedSeqs, ptddg::PITAOptions, TestIOPITA> TestSeedPITA;
typedef TestSite<ptddg::PITASeedSites, ptddg::PITASeedSeqs, ptddg::PITAOptions, TestIOPITA> TestSitePITA;

#endif //MIKAN_TEST_PITA_HPP_
