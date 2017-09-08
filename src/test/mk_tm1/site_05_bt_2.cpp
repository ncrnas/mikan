#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site05BT2 : public TestSiteTM1 {
protected:
    Site05BT2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_2.fasta";
    }

};

TEST_F(Site05BT2, mir1_bt) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(7u, sites.get_length());

    test_sites(sites, 0, "6mer", 5, 25, true, 0);
    test_sites(sites, 1, "7mer-m8", 10, 24, true, 0);
    test_sites(sites, 2, "7mer-m8", 24, 24, true, 0);
    test_sites(sites, 3, "6mer", 20, 24, true, 0);
    test_sites(sites, 4, "6mer", 21, 24, true, 0);
    test_sites(sites, 5, "6mer", 22, 24, true, 0);
    test_sites(sites, 6, "6mer", 23, 24, true, 0);


}

}
