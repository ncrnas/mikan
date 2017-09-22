#ifndef MIKAN_TEST_MIRANDA_HPP_
#define MIKAN_TEST_MIRANDA_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mk_main.hpp"
#include "mr3_core.hpp"
#include "mk_input.hpp"

typedef TestSeed<mr3as::MR3Options, mr3as::MR3SeedSeqs> TestSeedMR3AS;
typedef TestSite<mr3as::MR3Options, mr3as::MR3SeedSeqs, mr3as::MR3SeedSites> TestSiteMR3AS;

#endif //MIKAN_TEST_MIRANDA_HPP_
