#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site01Nmer1 : public TestSiteTM1 {
protected:
    Site01Nmer1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";
    }

};

TEST_F(Site01Nmer1, mir124) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(55u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 0, false, 0);
    test_sites(sites, 1, "6mer", 1, 1, false, 0);
    test_sites(sites, 2, "6mer", 2, 13, false, 0);
    test_sites(sites, 3, "6mer", 3, 14, false, 0);
    test_sites(sites, 4, "6mer", 4, 15, false, 0);
    test_sites(sites, 5, "6mer", 5, 18, false, 0);
    test_sites(sites, 6, "6mer", 6, 19, false, 0);
    test_sites(sites, 7, "6mer", 7, 20, true, 0);
    test_sites(sites, 8, "6mer", 8, 32, true, 0);
    test_sites(sites, 9, "6mer", 9, 33, true, 0);
    test_sites(sites, 10, "6mer", 10, 34, true, 0);

    test_sites(sites, 11, "7mer-A1", 11, 0, false, 0);
    test_sites(sites, 12, "7mer-A1", 12, 1, false, 0);
    test_sites(sites, 13, "7mer-A1", 13, 13, false, 0);
    test_sites(sites, 14, "7mer-A1", 14, 14, false, 0);
    test_sites(sites, 15, "7mer-A1", 15, 15, false, 0);
    test_sites(sites, 16, "7mer-A1", 16, 18, false, 0);
    test_sites(sites, 17, "7mer-A1", 17, 19, true, 0);
    test_sites(sites, 18, "7mer-A1", 18, 20, true, 0);
    test_sites(sites, 19, "7mer-A1", 19, 32, true, 0);
    test_sites(sites, 20, "7mer-A1", 20, 33, true, 0);
    test_sites(sites, 21, "6mer", 21, 34, true, 0);

    test_sites(sites, 22, "7mer-m8", 22, 1, false, 0);
    test_sites(sites, 23, "7mer-m8", 23, 13, false, 0);
    test_sites(sites, 24, "7mer-m8", 24, 14, false, 0);
    test_sites(sites, 25, "7mer-m8", 25, 15, false, 0);
    test_sites(sites, 26, "7mer-m8", 26, 18, false, 0);
    test_sites(sites, 27, "7mer-m8", 27, 19, false, 0);
    test_sites(sites, 28, "7mer-m8", 28, 20, true, 0);
    test_sites(sites, 29, "7mer-m8", 29, 21, true, 0);
    test_sites(sites, 30, "7mer-m8", 30, 32, true, 0);
    test_sites(sites, 31, "7mer-m8", 31, 33, true, 0);
    test_sites(sites, 32, "7mer-m8", 32, 34, true, 0);

    test_sites(sites, 33, "8mer", 33, 1, false, 0);
    test_sites(sites, 34, "8mer", 34, 13, false, 0);
    test_sites(sites, 35, "8mer", 35, 14, false, 0);
    test_sites(sites, 36, "8mer", 36, 15, false, 0);
    test_sites(sites, 37, "8mer", 37, 18, false, 0);
    test_sites(sites, 38, "8mer", 38, 19, true, 0);
    test_sites(sites, 39, "8mer", 39, 20, true, 0);
    test_sites(sites, 40, "8mer", 40, 21, true, 0);
    test_sites(sites, 41, "8mer", 41, 32, true, 0);
    test_sites(sites, 42, "8mer", 42, 33, true, 0);
    test_sites(sites, 43, "7mer-m8", 43, 34, true, 0);

    test_sites(sites, 44, "7mer-m8", 44, 2, false, 0);
    test_sites(sites, 45, "7mer-m8", 45, 13, false, 0);
    test_sites(sites, 46, "7mer-m8", 46, 14, false, 0);
    test_sites(sites, 47, "7mer-m8", 47, 15, false, 0);
    test_sites(sites, 48, "7mer-m8", 48, 18, false, 0);
    test_sites(sites, 49, "7mer-m8", 49, 19, false, 0);
    test_sites(sites, 50, "7mer-m8", 50, 20, true, 0);
    test_sites(sites, 51, "7mer-m8", 51, 21, true, 0);
    test_sites(sites, 52, "7mer-m8", 52, 32, true, 0);
    test_sites(sites, 53, "7mer-m8", 53, 33, true, 0);
    test_sites(sites, 54, "7mer-m8", 54, 34, true, 0);
}

}
