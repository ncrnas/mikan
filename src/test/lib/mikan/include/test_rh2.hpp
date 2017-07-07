#ifndef MIKAN_TEST_RH2_HPP_
#define MIKAN_TEST_RH2_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "rh2_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput> TestIORH2;
typedef TestSeed<rh2mfe::RH2SeedSeqs, TestIORH2> TestSeedRH2;
typedef TestSite<rh2mfe::RH2SeedSites, TestIORH2> TestSiteRH2;

#endif //MIKAN_TEST_RH2_HPP_
