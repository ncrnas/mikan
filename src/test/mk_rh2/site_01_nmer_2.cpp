#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site01Nmer2 : public TestSiteRH2 {
protected:
    Site01Nmer2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_2.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mer";
        mOverlapDef = "orig";
    }

};

TEST_F(Site01Nmer2, mir1_6mer) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(5u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 25, true, 0);
    test_sites(sites, 1, "6mer", 1, 25, true, 0);
    test_sites(sites, 2, "6mer", 2, 25, true, 0);
    test_sites(sites, 3, "6mer", 3, 25, true, 0);
    test_sites(sites, 4, "6mer", 4, 25, true, 0);
}

TEST_F(Site01Nmer2, mir1_7mer) {
    mSeedDef[0] = "7mer";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(3u, sites.get_length());

//    test_sites(sites, 0, "7mer", 0, 25, false);
//    test_sites(sites, 1, "7mer", 1, 25, false);
    test_sites(sites, 0, "7mer", 2, 25, true, 0);
    test_sites(sites, 1, "7mer", 3, 25, true, 0);
    test_sites(sites, 2, "7mer", 4, 25, true, 0);
}

TEST_F(Site01Nmer2, mir1_def) {
    mSeedDef[0] = "7mGU+";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(3u, sites.get_length());

//    test_sites(sites, 0, "7mer", 0, 25, false);
//    test_sites(sites, 1, "7mer", 1, 25, false);
    test_sites(sites, 0, "7mer", 2, 25, true, 0);
    test_sites(sites, 1, "7mer", 3, 25, true, 0);
    test_sites(sites, 2, "7mer", 4, 25, true, 0);
}
}
