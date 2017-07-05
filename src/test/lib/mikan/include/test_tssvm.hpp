#ifndef MIKAN_TEST_TSSVM_HPP_
#define MIKAN_TEST_TSSVM_HPP_

#include "test_main_io.hpp"
#include "test_seed.hpp"
#include "test_site.hpp"
#include "tssvm_core.hpp"

class TestIOTSSVM : public TestIOBase<tssvm::TSSVMCoreInput<tssvm::TRNATYPE>, tssvm::TSSVMOptions> {
    virtual void SetUp() {
        TestIOBase<tssvm::TSSVMCoreInput<tssvm::TRNATYPE>, tssvm::TSSVMOptions>::SetUp();

        modpath = STRINGIZE(TSSVM_MODEL_PATH);

        argc = 5;
        argv[1] = seqan::toCString(ifile1);
        argv[2] = seqan::toCString(ifile2);
        argv[3] = seqan::toCString(o2file1);
        argv[4] = seqan::toCString(o2file2);

    }

    seqan::CharString modpath;
};

typedef TestSeed<tssvm::TSSVMSeedSeqs<tssvm::TRNATYPE>, TestIOTSSVM> TestSeedTSSVM;
typedef TestSite<tssvm::TSSVMSeedSites<tssvm::TRNATYPE>, TestIOTSSVM> TestSiteTSSVM;

#endif //MIKAN_TEST_TSSVM_HPP_
