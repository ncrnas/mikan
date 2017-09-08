#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site04MMGU2 : public TestSiteTSSVM {
protected:
    Site04MMGU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_2.fasta";
    }

};

TEST_F(Site04MMGU2, mir1_mmgu) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(18u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);
    test_sites(sites, 2, "GUM", 23, 24, true, 1);
    test_sites(sites, 3, "GUM", 17, 24, true, 2);
    test_sites(sites, 4, "GUT", 11, 24, true, 5);
    test_sites(sites, 5, "GUM", 5, 24, true, 6);
    test_sites(sites, 6, "BT", 7, 23, true, 6);
    test_sites(sites, 7, "BT", 82, 24, true, 5);
    test_sites(sites, 8, "BT", 84, 24, true, 5);
    test_sites(sites, 9, "BT", 16, 24, true, 5);
    test_sites(sites, 10, "BT", 86, 24, true, 4);
    test_sites(sites, 11, "BT", 88, 24, true, 4);
    test_sites(sites, 12, "BT", 16, 24, true, 4);
    test_sites(sites, 13, "BT", 16, 24, true, 3);
    test_sites(sites, 14, "BT", 22, 24, true, 2);
    test_sites(sites, 15, "BM", 27, 24, true, 6);
    test_sites(sites, 16, "BM", 82, 25, true, 5);
    test_sites(sites, 17, "BM", 84, 25, true, 5);
}

}
