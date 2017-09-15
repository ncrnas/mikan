#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class Site01Nmer2 : public TestSiteTS5 {
protected:
    Site01Nmer2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_2.fasta";
    }

};

TEST_F(Site01Nmer2, mir1) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(4u, sites.get_length());

    test_sites(sites, 0, "7mer-A1", 1, 25, true, 0);
    test_sites(sites, 1, "7mer-m8", 2, 25, true, 0);
    test_sites(sites, 2, "8mer", 3, 25, true, 0);
    test_sites(sites, 3, "7mer-m8", 4, 25, true, 0);
}

}
