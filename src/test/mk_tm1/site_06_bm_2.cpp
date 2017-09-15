#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site06BM2 : public TestSiteTM1 {
protected:
    Site06BM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_2.fasta";
    }

};

TEST_F(Site06BM2, mir1_bm) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(2u, sites.get_length());

    test_sites(sites, 0, "6mer", 7, 26, true, 0);
    test_sites(sites, 1, "6mer", 6, 26, true, 0);
}

}
