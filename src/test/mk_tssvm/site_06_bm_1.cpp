#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

class Site06BM1 : public TestSiteTSSVM {
protected:
    Site06BM1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_1.fasta";
    }

};

TEST_F(Site06BM1, mir124_bm) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(8u, sites.get_length());

    test_sites(sites, 0, "BM", 1, 25, true, 6);
    test_sites(sites, 1, "BM", 3, 25, true, 5);
    test_sites(sites, 2, "BM", 5, 25, true, 4);
    test_sites(sites, 3, "BM", 5, 25, true, 3);
    test_sites(sites, 4, "BM", 7, 25, true, 2);
    test_sites(sites, 5, "BM", 7, 25, true, 1);

    test_sites(sites, 6, "LP", 7, 26, true, 1);
    test_sites(sites, 7, "LP", 6, 26, true, 1);
}

}
