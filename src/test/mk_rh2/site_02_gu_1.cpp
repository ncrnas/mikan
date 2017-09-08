#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site02GU1 : public TestSiteRH2 {
protected:
    Site02GU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_1.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mGU1";
        mOverlapDef = "orig";
    }

};

TEST_F(Site02GU1, mir124_6mer_gu1) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(11u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);

    test_sites(sites, 1, "6mer_GUT", 1, 24, true, 0);
    test_sites(sites, 2, "6mer_GUT", 2, 24, true, 0);
    test_sites(sites, 3, "6mer_GUT", 3, 24, true, 0);
    test_sites(sites, 4, "6mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "6mer_GUT", 5, 24, true, 0);

    test_sites(sites, 6, "6mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "6mer_GUT", 7, 24, true, 0);
    test_sites(sites, 8, "6mer_GUT", 8, 24, true, 0);
    test_sites(sites, 9, "6mer_GUT", 9, 24, true, 0);
    test_sites(sites, 10, "6mer_GUT", 10, 24, true, 0);
}

TEST_F(Site02GU1, mir124_6mer_gu_plus) {
    mSeedDef[0] = "6mGU+";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(16u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);

    test_sites(sites, 1, "6mer_GUT", 1, 24, true, 0);
    test_sites(sites, 2, "6mer_GUT", 2, 24, true, 0);
    test_sites(sites, 3, "6mer_GUT", 3, 24, true, 0);
    test_sites(sites, 4, "6mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "6mer_GUT", 5, 24, true, 0);

    test_sites(sites, 6, "6mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "6mer_GUT", 7, 24, true, 0);
    test_sites(sites, 8, "6mer_GUT", 8, 24, true, 0);
    test_sites(sites, 9, "6mer_GUT", 9, 24, true, 0);
    test_sites(sites, 10, "6mer_GUT", 10, 24, true, 0);

    test_sites(sites, 11, "6mer_GUT", 11, 24, true, 0);
    test_sites(sites, 12, "6mer_GUT", 12, 24, true, 0);
    test_sites(sites, 13, "6mer_GUT", 13, 24, true, 0);
    test_sites(sites, 14, "6mer_GUT", 14, 24, true, 0);
    test_sites(sites, 15, "6mer_GUT", 15, 24, true, 0);
}

TEST_F(Site02GU1, mir124_7mer_gu1) {
    mSeedDef[0] = "7mGU1";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(9u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);
//    test_sites(sites, 1, "7mer_GUT", 1, 24, false);
    test_sites(sites, 1, "7mer_GUT", 2, 24, true, 0);
    test_sites(sites, 2, "7mer_GUT", 3, 24, true, 0);
    test_sites(sites, 3, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 4, "7mer_GUT", 5, 24, true, 0);
//    test_sites(sites, 6, "7mer_GUT", 6, 24, false);
    test_sites(sites, 5, "7mer_GUT", 7, 24, true, 0);
    test_sites(sites, 6, "7mer_GUT", 8, 24, true, 0);
    test_sites(sites, 7, "7mer_GUT", 9, 24, true, 0);
    test_sites(sites, 8, "7mer_GUT", 10, 24, true, 0);
}

TEST_F(Site02GU1, mir124_def) {
    mSeedDef[0] = "7mGU+";
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);
    
    EXPECT_EQ(13u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);
//    test_sites(sites, 1, "7mer_GUT", 1, 24, false);
    test_sites(sites, 1, "7mer_GUT", 2, 24, true, 0);
    test_sites(sites, 2, "7mer_GUT", 3, 24, true, 0);
    test_sites(sites, 3, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 4, "7mer_GUT", 5, 24, true, 0);
//    test_sites(sites, 6, "7mer_GUT", 6, 24, false);
    test_sites(sites, 5, "7mer_GUT", 7, 24, true, 0);
    test_sites(sites, 6, "7mer_GUT", 8, 24, true, 0);
    test_sites(sites, 7, "7mer_GUT", 9, 24, true, 0);
    test_sites(sites, 8, "7mer_GUT", 10, 24, true, 0);
//    test_sites(sites, 11, "7mer_GUT", 11, 24, false);
    test_sites(sites, 9, "7mer_GUT", 12, 24, true, 0);
    test_sites(sites, 10, "7mer_GUT", 13, 24, true, 0);
    test_sites(sites, 11, "7mer_GUT", 14, 24, true, 0);
    test_sites(sites, 12, "7mer_GUT", 15, 24, true, 0);
}
}
