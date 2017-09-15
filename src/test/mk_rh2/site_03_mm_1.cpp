#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site03MM1 : public TestSiteRH2 {
protected:
    Site03MM1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_03_mm_1.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mGU1";
        mOverlapDef = "orig";
    }

};

TEST_F(Site03MM1, mir124_mm) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(4u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "6mer", 2, 24, true, 0);
    test_sites(sites, 3, "6mer", 3, 24, true, 0);
}

TEST_F(Site03MM1, mir124_def) {
    mSeedDef[0] = "7mGU+";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(0u, sites.get_length());

//    test_sites(sites, 0, "7mer", 0, 24, false);
//    test_sites(sites, 1, "7mer", 1, 24, false);
//    test_sites(sites, 2, "7mer", 2, 24, false);
//    test_sites(sites, 3, "7mer", 3, 24, false);
}
}
