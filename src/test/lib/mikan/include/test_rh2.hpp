#ifndef MIKAN_TEST_RH2_HPP_
#define MIKAN_TEST_RH2_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "rh2_core.hpp"

typedef TestIOBase<rh2mfe::RH2CoreInput<rh2mfe::TRNATYPE>, rh2mfe::RH2Options> TestIORH2;
typedef TestSeed<rh2mfe::RH2SeedSeqs<rh2mfe::TRNATYPE>, TestIORH2> TestSeedRH2;
typedef TestSite<rh2mfe::RH2SeedSites<rh2mfe::TRNATYPE>, TestIORH2> TestSiteRH2;

#endif //MIKAN_TEST_RH2_HPP_
