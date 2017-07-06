#ifndef MIKAN_TEST_TM1_HPP_
#define MIKAN_TEST_TM1_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "tm1_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput<tm1p::TRNATYPE> > TestIOTM1;
typedef TestSeed<tm1p::TM1SeedSeqs<tm1p::TRNATYPE>, TestIOTM1> TestSeedTM1;
typedef TestSite<tm1p::TM1SeedSites<tm1p::TRNATYPE>, TestIOTM1> TestSiteTM1;

#endif //MIKAN_TEST_TM1_HPP_
