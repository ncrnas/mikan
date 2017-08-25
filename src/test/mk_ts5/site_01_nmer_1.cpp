#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class Site01Nmer1 : public TestSiteTS5 {
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
    
    EXPECT_EQ(31u, sites.get_length());

    test_sites(sites, 1, "7mer-A1", 15, 15, true, 0);
    test_sites(sites, 2, "7mer-A1", 16, 18, true, 0);
    test_sites(sites, 3, "7mer-A1", 17, 19, true, 0);
    test_sites(sites, 4, "7mer-A1", 18, 20, true, 0);
    test_sites(sites, 5, "7mer-A1", 19, 32, true, 0);
    test_sites(sites, 6, "7mer-A1", 20, 33, true, 0);
    test_sites(sites, 7, "7mer-m8", 25, 15, true, 0);
    test_sites(sites, 8, "7mer-m8", 26, 18, true, 0);
    test_sites(sites, 9, "7mer-m8", 27, 19, true, 0);
    test_sites(sites, 10, "7mer-m8", 28, 20, true, 0);
    test_sites(sites, 11, "7mer-m8", 29, 21, true, 0);
    test_sites(sites, 12, "7mer-m8", 30, 32, true, 0);
    test_sites(sites, 13, "7mer-m8", 31, 33, true, 0);
    test_sites(sites, 14, "7mer-m8", 32, 34, true, 0);
    test_sites(sites, 15, "8mer", 36, 15, true, 0);
    test_sites(sites, 16, "8mer", 37, 18, true, 0);
    test_sites(sites, 17, "8mer", 38, 19, true, 0);
    test_sites(sites, 18, "8mer", 39, 20, true, 0);
    test_sites(sites, 19, "8mer", 40, 21, true, 0);
    test_sites(sites, 20, "8mer", 41, 32, true, 0);
    test_sites(sites, 21, "8mer", 42, 33, true, 0);
    test_sites(sites, 22, "7mer-m8", 43, 34, true, 0);
    test_sites(sites, 23, "7mer-m8", 47, 15, true, 0);
    test_sites(sites, 24, "7mer-m8", 48, 18, true, 0);
    test_sites(sites, 25, "7mer-m8", 49, 19, true, 0);
    test_sites(sites, 26, "7mer-m8", 50, 20, true, 0);
    test_sites(sites, 27, "7mer-m8", 51, 21, true, 0);
    test_sites(sites, 28, "7mer-m8", 52, 32, true, 0);
    test_sites(sites, 29, "7mer-m8", 53, 33, true, 0);
    test_sites(sites, 30, "7mer-m8", 54, 34, true, 0);
}

}
