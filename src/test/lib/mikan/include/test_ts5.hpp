#ifndef MIKAN_TEST_TS5_HPP_
#define MIKAN_TEST_TS5_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "ts5_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput<mikan::TRNATYPE> > TestIOTS5;
typedef TestSeed<ts5cs::TS5SeedSeqs<mikan::TRNATYPE>, TestIOTS5> TestSeedTS5;
typedef TestSite<ts5cs::TS5SeedSites<mikan::TRNATYPE>, TestIOTS5> TestSiteTS5;

#endif //MIKAN_TEST_TS5_HPP_
