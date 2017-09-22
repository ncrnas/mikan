#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site04MMGU1 : public TestSiteTSSVM {
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

    EXPECT_EQ(10u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);
    test_sites(sites, 2, "GUM", 11, 24, true, 3);
    test_sites(sites, 3, "GUM", 5, 24, true, 4);
    test_sites(sites, 4, "LP", 25, 24, true, 1);
    test_sites(sites, 5, "LP", 23, 24, true, 2);
    test_sites(sites, 6, "LP", 21, 24, true, 3);
    test_sites(sites, 7, "LP", 19, 24, true, 4);
    test_sites(sites, 8, "LP", 17, 24, true, 5);
    test_sites(sites, 9, "LP", 15, 24, true, 6);

}

}
