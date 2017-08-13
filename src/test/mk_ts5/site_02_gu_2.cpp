#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class Site02GU2 : public TestSiteTS5 {
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
    
    EXPECT_EQ(2u, sites.get_length());

    test_sites(sites, 0, "7mer-A1", 1, 24, true, 0);
    test_sites(sites, 1, "7mer-A1", 3, 24, true, 0);

}
}
