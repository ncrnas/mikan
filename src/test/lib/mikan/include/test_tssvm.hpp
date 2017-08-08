#ifndef MIKAN_TEST_TSSVM_HPP_
#define MIKAN_TEST_TSSVM_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "tssvm_core.hpp"
#include "mk_input.hpp"

typedef TestIOBase<mikan::MKInput> TestIOTSSVM;
typedef TestSeed<tssvm::TSSVMSeedSeqs, tssvm::TSSVMOptions, TestIOTSSVM> TestSeedTSSVM;
typedef TestSite<tssvm::TSSVMSeedSites, tssvm::TSSVMSeedSeqs, tssvm::TSSVMOptions, TestIOTSSVM> TestSiteTSSVM;

#endif //MIKAN_TEST_TSSVM_HPP_
