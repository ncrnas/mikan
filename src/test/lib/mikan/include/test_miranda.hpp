#ifndef MIKAN_TEST_MIRANDA_HPP_
#define MIKAN_TEST_MIRANDA_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mr3_core.hpp"

typedef TestIOBase<mr3as::MR3CoreInput<mr3as::TRNATYPE>, mr3as::MR3Options> TestIOMR3AS;
typedef TestSeed<mr3as::MR3SeedSeqs<mr3as::TRNATYPE>, TestIOMR3AS> TestSeedMR3AS;
typedef TestSite<mr3as::MR3SeedSites<mr3as::TRNATYPE>, TestIOMR3AS> TestSiteMR3AS;

#endif //MIKAN_TEST_MIRANDA_HPP_