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
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        mSeedDef1 = "6mer";
        mOverlapDef = "orig";
    }

    typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TIndexQGram TIdx;
    typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TFinder TFin;
    typedef rh2mfe::RH2SeedSites<rh2mfe::TRNATYPE> TSit;

};

TEST_F(Site05BT2, mir1_bt) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef1, mOverlapDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(0u, sites.get_length());
}

TEST_F(Site05BT2, mir1_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef1 = "7mGU+";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef1, mOverlapDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(28u, sites.get_length());

    test_sites2(sites, 0, "7mer_GUT", 5, 25, false);

    test_sites2(sites, 1, "7mer_GUT", 10, 24, true);

    test_sites2(sites, 2, "7mer_GUT", 24, 24, true);

    test_sites2(sites, 3, "7mer_GUT", 0, 15, true);
    test_sites2(sites, 4, "7mer_GUT", 1, 15, true);
    test_sites2(sites, 5, "7mer_GUT", 2, 15, true);
    test_sites2(sites, 6, "7mer_GUT", 3, 15, true);
    test_sites2(sites, 7, "7mer_GUT", 4, 15, true);
    test_sites2(sites, 8, "7mer_GUT", 5, 15, true);
    test_sites2(sites, 9, "7mer_GUT", 6, 15, true);
    test_sites2(sites, 10, "7mer_GUT", 7, 15, true);
    test_sites2(sites, 11, "7mer_GUT", 8, 15, true);
    test_sites2(sites, 12, "7mer_GUT", 9, 15, true);
    test_sites2(sites, 13, "7mer_GUT", 10, 15, true);
    test_sites2(sites, 14, "7mer_GUT", 11, 15, true);
    test_sites2(sites, 15, "7mer_GUT", 12, 15, true);
    test_sites2(sites, 16, "7mer_GUT", 13, 15, true);
    test_sites2(sites, 17, "7mer_GUT", 14, 15, true);
    test_sites2(sites, 18, "7mer_GUT", 15, 15, true);
    test_sites2(sites, 19, "7mer_GUT", 16, 15, true);
    test_sites2(sites, 20, "7mer_GUT", 17, 15, true);
    test_sites2(sites, 21, "7mer_GUT", 18, 15, true);
    test_sites2(sites, 22, "7mer_GUT", 19, 15, true);
    test_sites2(sites, 23, "7mer_GUT", 20, 15, true);
    test_sites2(sites, 24, "7mer_GUT", 21, 15, true);
    test_sites2(sites, 25, "7mer_GUT", 22, 15, true);
    test_sites2(sites, 26, "7mer_GUT", 23, 15, true);
    test_sites2(sites, 27, "7mer_GUT", 24, 15, true);
}

}
