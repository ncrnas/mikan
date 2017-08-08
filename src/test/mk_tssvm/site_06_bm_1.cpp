#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site06BM1 : public TestSiteTSSVM {
protected:
    Site06BM1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tssvm/";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef tssvm::TSSVMSeedSites TSit;
    typedef tssvm::TSSVMSeedSeqs TSeed;

};

TEST_F(Site06BM1, mir124_bm) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;


    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(8u, sites.get_length());

    test_sites(sites, 0, "BM", 1, 25, true, 6);
    test_sites(sites, 1, "BM", 3, 25, true, 5);
    test_sites(sites, 2, "BM", 5, 25, true, 4);
    test_sites(sites, 3, "BM", 5, 25, true, 3);
    test_sites(sites, 4, "BM", 7, 25, true, 2);
    test_sites(sites, 5, "BM", 7, 25, true, 1);

    test_sites(sites, 6, "LP", 7, 26, true, 1);
    test_sites(sites, 7, "LP", 6, 26, true, 1);
}

}
