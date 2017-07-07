#ifndef MIKAN_TEST_TSSVM_HPP_
#define MIKAN_TEST_TSSVM_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "tssvm_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput<mikan::TRNATYPE> > TestIOTSSVM;
typedef TestSeed<tssvm::TSSVMSeedSeqs<mikan::TRNATYPE>, TestIOTSSVM> TestSeedTSSVM;
typedef TestSite<tssvm::TSSVMSeedSites<mikan::TRNATYPE>, TestIOTSSVM> TestSiteTSSVM;

#endif //MIKAN_TEST_TSSVM_HPP_
