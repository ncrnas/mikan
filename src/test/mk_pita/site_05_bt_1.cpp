#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site05BT1 : public TestSitePITA {
protected:
    Site05BT1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_1.fasta";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "0";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "1";
    }

};

TEST_F(Site05BT1, mir124_bt) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(0u, sites.get_length());
}

TEST_F(Site05BT1, mir124_def) {
    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(10u, sites.get_length());

//    test_sites(sites, 0, "MM", 4, 25, false, 0);
//    test_sites(sites, 1, "MM", 5, 25, false, 0);

//    test_sites(sites, 2, "MM", 0, 25, false, 0);
//    test_sites(sites, 3, "MM", 1, 25, false, 0);
//    test_sites(sites, 4, "MM", 2, 25, false, 0);
//    test_sites(sites, 5, "MM", 3, 25, false, 0);

    test_sites(sites, 0, "8mer_MM", 9, 24, true, 4);

    test_sites(sites, 1, "8mer_MM", 19, 24, true, 4);

//    test_sites(sites, 8, "MM", 16, 24, false, 0);

    test_sites(sites, 2, "8mer_MM", 17, 24, true, 4);
    test_sites(sites, 3, "8mer_MM", 18, 24, true, 4);

    test_sites(sites, 4, "8mer_MM", 23, 24, true, 5);
    test_sites(sites, 5, "8mer_MM", 24, 24, true, 5);

//    test_sites(sites, 13, "MM", 20, 24, false, 0);
    test_sites(sites, 6, "8mer_MM", 21, 24, true, 5);
    test_sites(sites, 7, "8mer_MM", 22, 24, true, 5);

    test_sites(sites, 8, "8mer_MMGU", 10, 24, true, 4);

    test_sites(sites, 9, "8mer_MMGU", 15, 24, true, 4);
}

}
