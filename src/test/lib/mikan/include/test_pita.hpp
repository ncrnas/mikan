#ifndef MIKAN_TEST_PITA_HPP_
#define MIKAN_TEST_PITA_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "pita_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput<mikan::TRNATYPE> > TestIOPITA;
typedef TestSeed<ptddg::PITASeedSeqs<mikan::TRNATYPE>, TestIOPITA> TestSeedPITA;
typedef TestSite<ptddg::PITASeedSites<mikan::TRNATYPE>, TestIOPITA> TestSitePITA;

#endif //MIKAN_TEST_PITA_HPP_
