#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site03MM2 : public TestSiteTSSVM {
protected:
    Site03MM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_03_mm_2.fasta";
    }

};

TEST_F(Site03MM2, mir1_mm) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(35u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);
    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "7mer-A1", 3, 24, true, 0);

//        test_sites(sites, 4, "BT", 1, 23, true, 7);
//        test_sites(sites, 5, "BT", 5, 23, true, 6);
//        test_sites(sites, 6, "BT", 1, 23, true, 6);
//        test_sites(sites, 7, "BT", 5, 23, true, 5);
//        test_sites(sites, 8, "BT", 30, 24, true, 3);
//        test_sites(sites, 9, "BT", 32, 24, true, 3);
//        test_sites(sites, 10, "BT", 36, 24, true, 2);
//        test_sites(sites, 11, "BT", 38, 24, true, 2);
//
//        test_sites(sites, 12, "BM", 5, 24, true, 6);
//        test_sites(sites, 13, "BM", 7, 24, true, 6);
//        test_sites(sites, 14, "BM", 9, 24, true, 6);
//        test_sites(sites, 15, "BM", 36, 23, true, 2);
//        test_sites(sites, 16, "BM", 37, 23, true, 2);
//        test_sites(sites, 17, "BM", 38, 23, true, 2);
//        test_sites(sites, 18, "BM", 39, 23, true, 2);
//        test_sites(sites, 19, "BM", 36, 23, true, 1);
//        test_sites(sites, 20, "BM", 37, 23, true, 1);
//        test_sites(sites, 21, "BM", 38, 23, true, 1);
//        test_sites(sites, 22, "BM", 39, 23, true, 1);

    test_sites(sites, 23, "LP", 37, 24, true, 1);
    test_sites(sites, 24, "LP", 39, 24, true, 1);
    test_sites(sites, 25, "LP", 31, 24, true, 2);
    test_sites(sites, 26, "LP", 33, 24, true, 2);
    test_sites(sites, 27, "LP", 25, 24, true, 3);
    test_sites(sites, 28, "LP", 27, 24, true, 3);
    test_sites(sites, 29, "LP", 19, 24, true, 4);
    test_sites(sites, 30, "LP", 21, 24, true, 4);
    test_sites(sites, 31, "LP", 13, 24, true, 5);
    test_sites(sites, 32, "LP", 15, 24, true, 5);
    test_sites(sites, 33, "LP", 7, 24, true, 6);
    test_sites(sites, 34, "LP", 9, 24, true, 6);

}

}
