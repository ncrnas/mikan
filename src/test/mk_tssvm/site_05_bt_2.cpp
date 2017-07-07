#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site05BT2 : public TestSiteTSSVM {
protected:
    Site05BT2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tssvm/";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef tssvm::TSSVMSeedSites<mikan::TRNATYPE> TSit;

};

TEST_F(Site05BT2, mir1_bt) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[1]);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(26u, sites.get_length());

    test_sites(sites, 0, "BT", 1, 24, true, 6);
    test_sites(sites, 1, "BT", 3, 24, true, 6);

    test_sites(sites, 2, "BT", 1, 24, true, 5);
    test_sites(sites, 3, "BT", 3, 24, true, 5);
    test_sites(sites, 4, "BT", 8, 24, true, 5);
    test_sites(sites, 5, "BT", 9, 24, true, 5);
    test_sites(sites, 6, "BT", 10, 24, true, 5);

    test_sites(sites, 7, "BT", 13, 24, true, 4);
    test_sites(sites, 8, "BT", 14, 24, true, 4);
    test_sites(sites, 9, "BT", 15, 24, true, 4);
    test_sites(sites, 10, "BT", 10, 24, true, 4);

    test_sites(sites, 11, "BT", 18, 24, true, 3);
    test_sites(sites, 12, "BT", 19, 24, true, 3);
    test_sites(sites, 13, "BT", 10, 24, true, 3);

    test_sites(sites, 14, "BT", 22, 24, true, 2);
    test_sites(sites, 15, "BT", 23, 24, true, 2);
    test_sites(sites, 16, "BT", 24, 24, true, 2);

    test_sites(sites, 17, "BM", 1, 25, true, 6);
    test_sites(sites, 18, "BM", 3, 25, true, 6);
    test_sites(sites, 19, "BM", 8, 25, true, 5);
    test_sites(sites, 20, "BM", 20, 23, true, 2);
    test_sites(sites, 21, "BM", 21, 23, true, 2);
    test_sites(sites, 22, "BM", 22, 23, true, 2);
    test_sites(sites, 23, "BM", 20, 23, true, 1);
    test_sites(sites, 24, "BM", 21, 23, true, 1);
    test_sites(sites, 25, "BM", 22, 23, true, 1);
}

}
