#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site04MMGU1 : public TestSiteTM1 {
protected:
    Site04MMGU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_1.fasta";
    }

};

TEST_F(Site04MMGU1, mir124_mmgu) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(24u, sites.get_length());

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

    test_sites(sites, 14, "6mer", 24, 24, true, 0);
    test_sites(sites, 15, "6mer", 25, 24, true, 0);

    test_sites(sites, 16, "6mer", 42, 24, true, 0);
    test_sites(sites, 17, "6mer", 43, 24, true, 0);
    test_sites(sites, 18, "6mer", 44, 24, true, 0);
    test_sites(sites, 19, "6mer", 45, 24, true, 0);

    test_sites(sites, 20, "6mer", 62, 24, true, 0);
    test_sites(sites, 21, "6mer", 63, 24, true, 0);
    test_sites(sites, 22, "6mer", 64, 24, true, 0);
    test_sites(sites, 23, "6mer", 65, 24, true, 0);
}
}
