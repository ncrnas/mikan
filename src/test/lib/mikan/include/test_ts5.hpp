#ifndef MIKAN_TEST_TS5_HPP_
#define MIKAN_TEST_TS5_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "mk_core_main.hpp"
#include "ts5_core.hpp"
#include "mk_input.hpp"

typedef TestSeed<ts5cs::TS5Options, ts5cs::TS5SeedSeqs> TestSeedTS5;
typedef TestSite<ts5cs::TS5Options, ts5cs::TS5SeedSeqs, ts5cs::TS5SeedSites> TestSiteTS5;

#endif //MIKAN_TEST_TS5_HPP_
