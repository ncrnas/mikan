#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class Site06BM2 : public TestSiteMR3AS {
protected:
    Site06BM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_2.fasta";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "1";
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "1";
    }

};

TEST_F(Site06BM2, mir1_bm) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(2u, sites.get_length());

//    test_sites(sites, 0, "MM", 0, 25, false, 0);
//    test_sites(sites, 1, "MM", 1, 25, false, 0);
    test_sites(sites, 0, "8mer_MM", 7, 26, true, 5);
    test_sites(sites, 1, "8mer_MM", 6, 26, true, 5);

//    test_sites(sites, 4, "BT", 6, 25, false, 0);
//    test_sites(sites, 5, "BT", 7, 25, false, 0);
//    test_sites(sites, 6, "BT", 7, 26, false, 0);
//    test_sites(sites, 7, "BT", 6, 26, false, 0);
}

TEST_F(Site06BM2, mir1_def) {
    mSeedDef[3] = "+";
    mSeedDef[4] = "1:1";
    mSeedDef[5] = "1";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(10u, sites.get_length());

    test_sites(sites, 0, "7mer_GU+", 0, 15, true, 0);
    test_sites(sites, 1, "7mer_GU+", 1, 15, true, 0);
    test_sites(sites, 2, "7mer_GU+", 2, 15, true, 0);
    test_sites(sites, 3, "7mer_GU+", 3, 15, true, 0);
    test_sites(sites, 4, "7mer_GU+", 4, 15, true, 0);
    test_sites(sites, 5, "7mer_GU+", 5, 15, true, 0);
    test_sites(sites, 6, "7mer_GU+", 6, 15, true, 0);
    test_sites(sites, 7, "7mer_GU+", 7, 15, true, 0);
//    test_sites(sites, 8, "MM", 0, 25, false, 0);
//    test_sites(sites, 9, "MM", 1, 25, false, 0);
    test_sites(sites, 8, "8mer_MM", 7, 26, true, 5);
    test_sites(sites, 9, "8mer_MM", 6, 26, true, 5);

//    test_sites(sites, 12, "BT", 6, 25, false, 0);
//    test_sites(sites, 13, "BT", 7, 25, false, 0);
//    test_sites(sites, 14, "BT", 7, 26, false, 0);
//    test_sites(sites, 15, "BT", 6, 26, false, 0);
}
}
