#ifndef MIKAN_TEST_PITA_HPP_
#define MIKAN_TEST_PITA_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "pita_core.hpp"

typedef TestIOBase<ptddg::PITACoreInput<ptddg::TRNATYPE>, ptddg::PITAOptions> TestIOPITA;
typedef TestSeed<ptddg::PITASeedSeqs<ptddg::TRNATYPE>, TestIOPITA> TestSeedPITA;
typedef TestSite<ptddg::PITASeedSites<ptddg::TRNATYPE>, TestIOPITA> TestSitePITA;

#endif //MIKAN_TEST_PITA_HPP_
