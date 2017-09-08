#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site02GU1 : public TestSiteTSSVM {
protected:
    Site02GU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_1.fasta";
    }

};

TEST_F(Site02GU1, mir124_gu) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(5u, sites.get_length());

    test_sites(sites, 0, "7mer-m8", 0, 24, true, 0);

    test_sites(sites, 1, "GUM", 8, 24, true, 3);
    test_sites(sites, 2, "GUM", 10, 24, true, 3);
    test_sites(sites, 3, "GUM", 3, 24, true, 4);
    test_sites(sites, 4, "GUM", 5, 24, true, 4);
}
}
