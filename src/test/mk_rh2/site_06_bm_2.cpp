#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site06BM2 : public TestSiteRH2 {
protected:
    Site06BM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_2.fasta";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mer";
        mOverlapDef = "orig";
    }

};

TEST_F(Site06BM2, mir1_bm_6mer) {
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(0u, sites.get_length());
}

TEST_F(Site06BM2, mir1_def) {
    mSeedDef[0] = "7mGU+";
    create_seed_seqs(1);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(8u, sites.get_length());

    test_sites(sites, 0, "7mer_GUT", 0, 15, true, 0);
    test_sites(sites, 1, "7mer_GUT", 1, 15, true, 0);
    test_sites(sites, 2, "7mer_GUT", 2, 15, true, 0);
    test_sites(sites, 3, "7mer_GUT", 3, 15, true, 0);
    test_sites(sites, 4, "7mer_GUT", 4, 15, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 15, true, 0);
    test_sites(sites, 6, "7mer_GUT", 6, 15, true, 0);
    test_sites(sites, 7, "7mer_GUT", 7, 15, true, 0);
}

}
