#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site02GU1 : public TestSitePITA {
protected:
    Site02GU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_1.fasta";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "1";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "0";
    }

};

TEST_F(Site02GU1, mir124_gu) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(9u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);
//    test_sites(sites, 1, "GUT", 1, 24, false, 0);
    test_sites(sites, 1, "7mer_GUT", 2, 24, true, 2);
    test_sites(sites, 2, "7mer_GUT", 3, 24, true, 2);
    test_sites(sites, 3, "8mer_GUT", 4, 24, true, 2);
    test_sites(sites, 4, "8mer_GUT", 5, 24, true, 2);
//    test_sites(sites, 6, "GUT", 6, 24, false, 0);
    test_sites(sites, 5, "7mer_GUT", 7, 24, true, 3);
    test_sites(sites, 6, "7mer_GUT", 8, 24, true, 3);
    test_sites(sites, 7, "8mer_GUT", 9, 24, true, 3);
    test_sites(sites, 8, "8mer_GUT", 10, 24, true, 3);
}

TEST_F(Site02GU1, mir124_gu_plus) {
    mSeedDef[3] = "+";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(13u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);
//    test_sites(sites, 1, "GU+", 1, 24, false, 0);
    test_sites(sites, 1, "7mer_GU+", 2, 24, true, 0);
    test_sites(sites, 2, "7mer_GU+", 3, 24, true, 0);
    test_sites(sites, 3, "8mer_GU+", 4, 24, true, 0);
    test_sites(sites, 4, "8mer_GU+", 5, 24, true, 0);
//    test_sites(sites, 6, "GU+", 6, 24, false, 0);
    test_sites(sites, 5, "7mer_GU+", 7, 24, true, 0);
    test_sites(sites, 6, "7mer_GU+", 8, 24, true, 0);
    test_sites(sites, 7, "8mer_GU+", 9, 24, true, 0);
    test_sites(sites, 8, "8mer_GU+", 10, 24, true, 0);
//    test_sites(sites, 11, "GU+", 11, 24, false, 0);
    test_sites(sites, 9, "7mer_GU+", 12, 24, true, 0);
    test_sites(sites, 10, "7mer_GU+", 13, 24, true, 0);
    test_sites(sites, 11, "8mer_GU+", 14, 24, true, 0);
    test_sites(sites, 12, "8mer_GU+", 15, 24, true, 0);
}

TEST_F(Site02GU1, mir124_def) {
    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(9u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);

//    test_sites(sites, 1, "GUT", 1, 24, false, 0);
    test_sites(sites, 1, "7mer_GUT", 2, 24, true, 2);
    test_sites(sites, 2, "7mer_GUT", 3, 24, true, 2);
    test_sites(sites, 3, "8mer_GUT", 4, 24, true, 2);
    test_sites(sites, 4, "8mer_GUT", 5, 24, true, 2);

//    test_sites(sites, 5, "GUT", 6, 24, false, 0);
    test_sites(sites, 5, "7mer_GUT", 7, 24, true, 3);
    test_sites(sites, 6, "7mer_GUT", 8, 24, true, 3);
    test_sites(sites, 7, "8mer_GUT", 9, 24, true, 3);
    test_sites(sites, 8, "8mer_GUT", 10, 24, true, 3);

}
}
