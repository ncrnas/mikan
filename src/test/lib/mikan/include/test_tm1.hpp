#ifndef MIKAN_TEST_TM1_HPP_
#define MIKAN_TEST_TM1_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mk_main.hpp"
#include "tm1_core.hpp"
#include "mk_input.hpp"

typedef TestSeed<tm1p::TM1Options, tm1p::TM1SeedSeqs> TestSeedTM1;
typedef TestSite<tm1p::TM1Options, tm1p::TM1SeedSeqs, tm1p::TM1SeedSites> TestSiteTM1;

#endif //MIKAN_TEST_TM1_HPP_
