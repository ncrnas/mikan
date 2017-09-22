#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site04MMGU2 : public TestSiteTM1 {
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

    EXPECT_EQ(40u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);

    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "7mer-A1", 3, 24, true, 0);
    test_sites(sites, 4, "7mer-m8", 4, 24, true, 0);
    test_sites(sites, 5, "8mer", 5, 24, true, 0);
    test_sites(sites, 6, "6mer", 6, 24, true, 0);
    test_sites(sites, 7, "7mer-A1", 7, 24, true, 0);

    test_sites(sites, 8, "6mer", 8, 24, true, 0);
    test_sites(sites, 9, "7mer-A1", 9, 24, true, 0);
    test_sites(sites, 10, "7mer-m8", 10, 24, true, 0);
    test_sites(sites, 11, "8mer", 11, 24, true, 0);
    test_sites(sites, 12, "6mer", 12, 24, true, 0);
    test_sites(sites, 13, "7mer-A1", 13, 24, true, 0);

    test_sites(sites, 14, "6mer", 14, 24, true, 0);
    test_sites(sites, 15, "7mer-A1", 15, 24, true, 0);
    test_sites(sites, 16, "7mer-m8", 16, 24, true, 0);
    test_sites(sites, 17, "8mer", 17, 24, true, 0);
    test_sites(sites, 18, "6mer", 18, 24, true, 0);
    test_sites(sites, 19, "7mer-A1", 19, 24, true, 0);

    test_sites(sites, 20, "6mer", 20, 24, true, 0);
    test_sites(sites, 21, "7mer-A1", 21, 24, true, 0);
    test_sites(sites, 22, "7mer-m8", 22, 24, true, 0);
    test_sites(sites, 23, "8mer", 23, 24, true, 0);
    test_sites(sites, 24, "6mer", 24, 24, true, 0);
    test_sites(sites, 25, "7mer-A1", 25, 24, true, 0);

    test_sites(sites, 26, "6mer", 54, 24, true, 0);
    test_sites(sites, 27, "6mer", 55, 24, true, 0);
    test_sites(sites, 28, "6mer", 56, 24, true, 0);
    test_sites(sites, 29, "6mer", 57, 24, true, 0);

    test_sites(sites, 30, "6mer", 74, 24, true, 0);
    test_sites(sites, 31, "6mer", 75, 24, true, 0);
    test_sites(sites, 32, "6mer", 76, 24, true, 0);
    test_sites(sites, 33, "6mer", 77, 24, true, 0);

    test_sites(sites, 34, "6mer", 90, 24, true, 0);
    test_sites(sites, 35, "6mer", 91, 24, true, 0);
    test_sites(sites, 36, "6mer", 92, 24, true, 0);
    test_sites(sites, 37, "6mer", 93, 24, true, 0);

    test_sites(sites, 38, "6mer", 36, 24, true, 0);
    test_sites(sites, 39, "6mer", 37, 24, true, 0);
}

}
