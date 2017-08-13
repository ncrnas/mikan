#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site02GU2 : public TestSiteTM1 {
protected:
    Site02GU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_2.fasta";
    }

};

TEST_F(Site02GU2, mir1_gu) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(40u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);
    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "7mer-A1", 3, 24, true, 0);

    test_sites(sites, 4, "7mer-m8", 4, 24, true, 0);
    test_sites(sites, 5, "8mer", 5, 24, true, 0);
    test_sites(sites, 6, "7mer-m8", 6, 24, true, 0);
    test_sites(sites, 7, "8mer", 7, 24, true, 0);

    test_sites(sites, 8, "6mer", 8, 24, true, 0);
    test_sites(sites, 9, "6mer", 9, 24, true, 0);
    test_sites(sites, 10, "7mer-A1", 10, 24, true, 0);
    test_sites(sites, 11, "6mer", 11, 24, true, 0);
    test_sites(sites, 12, "7mer-A1", 12, 24, true, 0);

    test_sites(sites, 13, "7mer-m8", 13, 24, true, 0);
    test_sites(sites, 14, "8mer", 14, 24, true, 0);
    test_sites(sites, 15, "7mer-m8", 15, 24, true, 0);
    test_sites(sites, 16, "8mer", 16, 24, true, 0);

    test_sites(sites, 17, "6mer", 17, 24, true, 0);
    test_sites(sites, 18, "6mer", 18, 24, true, 0);
    test_sites(sites, 19, "7mer-A1", 19, 24, true, 0);
    test_sites(sites, 20, "6mer", 20, 24, true, 0);
    test_sites(sites, 21, "7mer-A1", 21, 24, true, 0);

    test_sites(sites, 22, "7mer-m8", 22, 24, true, 0);
    test_sites(sites, 23, "8mer", 23, 24, true, 0);
    test_sites(sites, 24, "7mer-m8", 24, 24, true, 0);
    test_sites(sites, 25, "8mer", 25, 24, true, 0);

    test_sites(sites, 26, "6mer", 26, 24, true, 0);
    test_sites(sites, 27, "6mer", 27, 24, true, 0);
    test_sites(sites, 28, "7mer-A1", 28, 24, true, 0);
    test_sites(sites, 29, "6mer", 29, 24, true, 0);
    test_sites(sites, 30, "7mer-A1", 30, 24, true, 0);

    test_sites(sites, 31, "7mer-m8", 31, 24, true, 0);
    test_sites(sites, 32, "8mer", 32, 24, true, 0);
    test_sites(sites, 33, "7mer-m8", 33, 24, true, 0);
    test_sites(sites, 34, "8mer", 34, 24, true, 0);

    test_sites(sites, 35, "6mer", 35, 24, true, 0);
    test_sites(sites, 36, "6mer", 36, 24, true, 0);
    test_sites(sites, 37, "7mer-A1", 37, 24, true, 0);
    test_sites(sites, 38, "6mer", 38, 24, true, 0);
    test_sites(sites, 39, "7mer-A1", 39, 24, true, 0);
}
}
