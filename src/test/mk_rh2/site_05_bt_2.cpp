#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site05BT2 : public TestSiteRH2 {
protected:
    Site05BT2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_2.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mer";
        mOverlapDef = "orig";
    }

};

TEST_F(Site05BT2, mir1_bt) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(0u, sites.get_length());
}

TEST_F(Site05BT2, mir1_def) {
    mSeedDef[0] = "7mGU+";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(27u, sites.get_length());

//    test_sites(sites, 0, "7mer_GUT", 5, 25, false);
    test_sites(sites, 0, "7mer_GUT", 10, 24, true, 0);
    test_sites(sites, 1, "7mer_GUT", 24, 24, true, 0);
    test_sites(sites, 2, "7mer_GUT", 0, 15, true, 0);
    test_sites(sites, 3, "7mer_GUT", 1, 15, true, 0);
    test_sites(sites, 4, "7mer_GUT", 2, 15, true, 0);
    test_sites(sites, 5, "7mer_GUT", 3, 15, true, 0);
    test_sites(sites, 6, "7mer_GUT", 4, 15, true, 0);
    test_sites(sites, 7, "7mer_GUT", 5, 15, true, 0);
    test_sites(sites, 8, "7mer_GUT", 6, 15, true, 0);
    test_sites(sites, 9, "7mer_GUT", 7, 15, true, 0);
    test_sites(sites, 10, "7mer_GUT", 8, 15, true, 0);
    test_sites(sites, 11, "7mer_GUT", 9, 15, true, 0);
    test_sites(sites, 12, "7mer_GUT", 10, 15, true, 0);
    test_sites(sites, 13, "7mer_GUT", 11, 15, true, 0);
    test_sites(sites, 14, "7mer_GUT", 12, 15, true, 0);
    test_sites(sites, 15, "7mer_GUT", 13, 15, true, 0);
    test_sites(sites, 16, "7mer_GUT", 14, 15, true, 0);
    test_sites(sites, 17, "7mer_GUT", 15, 15, true, 0);
    test_sites(sites, 18, "7mer_GUT", 16, 15, true, 0);
    test_sites(sites, 19, "7mer_GUT", 17, 15, true, 0);
    test_sites(sites, 20, "7mer_GUT", 18, 15, true, 0);
    test_sites(sites, 21, "7mer_GUT", 19, 15, true, 0);
    test_sites(sites, 22, "7mer_GUT", 20, 15, true, 0);
    test_sites(sites, 23, "7mer_GUT", 21, 15, true, 0);
    test_sites(sites, 24, "7mer_GUT", 22, 15, true, 0);
    test_sites(sites, 25, "7mer_GUT", 23, 15, true, 0);
    test_sites(sites, 26, "7mer_GUT", 24, 15, true, 0);

}

}
