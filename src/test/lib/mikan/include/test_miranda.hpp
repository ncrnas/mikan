#ifndef MIKAN_TEST_MIRANDA_HPP_
#define MIKAN_TEST_MIRANDA_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mr3_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput<mikan::TRNATYPE> > TestIOMR3AS;
typedef TestSeed<mr3as::MR3SeedSeqs<mikan::TRNATYPE>, TestIOMR3AS> TestSeedMR3AS;
typedef TestSite<mr3as::MR3SeedSites<mikan::TRNATYPE>, TestIOMR3AS> TestSiteMR3AS;

#endif //MIKAN_TEST_MIRANDA_HPP_
