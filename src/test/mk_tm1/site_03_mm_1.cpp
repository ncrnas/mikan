#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site03MM1 : public TestSiteTM1 {
protected:
    Site03MM1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_03_mm_1.fasta";
    }

};

TEST_F(Site03MM1, mir124_mm) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(8u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);
    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "7mer-A1", 3, 24, true, 0);

//    test_sites(sites, 4, "", 34, 24, false);
//    test_sites(sites, 5, "", 35, 24, false);
    test_sites(sites, 4, "6mer", 36, 24, true, 0);
    test_sites(sites, 5, "6mer", 37, 24, true, 0);
    test_sites(sites, 6, "6mer", 38, 24, true, 0);
    test_sites(sites, 7, "6mer", 39, 24, true, 0);

}

}
