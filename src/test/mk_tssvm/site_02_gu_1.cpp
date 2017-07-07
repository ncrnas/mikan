#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site02GU1 : public TestSiteTSSVM {
protected:
    Site02GU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tssvm/";
    }

    typedef tssvm::TSSVMCore<mikan::TRNATYPE>::TIndexQGram TIdx;
    typedef tssvm::TSSVMCore<mikan::TRNATYPE>::TFinder TFin;
    typedef tssvm::TSSVMSeedSites<mikan::TRNATYPE> TSit;

};

TEST_F(Site02GU1, mir124_gu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[0]);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(5u, sites.get_length());

    test_sites(sites, 0, "7mer-m8", 0, 24, true, 0);

    test_sites(sites, 1, "GUM", 8, 24, true, 3);
    test_sites(sites, 2, "GUM", 10, 24, true, 3);
    test_sites(sites, 3, "GUM", 3, 24, true, 4);
    test_sites(sites, 4, "GUM", 5, 24, true, 4);
}
}
