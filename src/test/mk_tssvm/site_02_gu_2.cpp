#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site02GU2 : public TestSiteTSSVM {
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

    EXPECT_EQ(21u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);
    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "7mer-A1", 3, 24, true, 0);

    test_sites(sites, 4, "GUM", 32, 24, true, 1);
    test_sites(sites, 5, "GUM", 34, 24, true, 1);
    test_sites(sites, 6, "GUM", 23, 24, true, 2);
    test_sites(sites, 7, "GUM", 25, 24, true, 2);
    test_sites(sites, 8, "GUT", 14, 24, true, 5);
    test_sites(sites, 9, "GUT", 16, 24, true, 5);
    test_sites(sites, 10, "GUM", 5, 24, true, 6);
    test_sites(sites, 11, "GUM", 7, 24, true, 6);
}

}
