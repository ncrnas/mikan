#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class Site05BT1 : public TestSiteMR3AS {
protected:
    Site05BT1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_miranda/";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "0";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "1";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef mr3as::MR3SeedSites TSit;
    typedef mr3as::MR3SeedSeqs TSeed;

};

TEST_F(Site05BT1, mir124_bt) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    TSeed seedSeqs;

    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(25u, sites.get_length());

//    test_sites(sites, 0, "BT", 20, 23, false, 0);
//    test_sites(sites, 1, "BT", 21, 23, false, 0);
//    test_sites(sites, 2, "BT", 22, 23, false, 0);
//    test_sites(sites, 3, "BT", 23, 23, false, 0);
//    test_sites(sites, 4, "BT", 24, 23, false, 0);

    test_sites(sites, 0, "7mer_BT", 0, 24, true, 1);
    test_sites(sites, 1, "7mer_BT", 1, 24, true, 1);
    test_sites(sites, 2, "8mer_BT", 2, 24, true, 1);
    test_sites(sites, 3, "8mer_BT", 3, 24, true, 1);
    test_sites(sites, 4, "7mer_BT", 4, 24, true, 1);
    test_sites(sites, 5, "7mer_BT", 5, 24, true, 1);

    test_sites(sites, 6, "7mer_BT", 6, 24, true, 2);
    test_sites(sites, 7, "8mer_BT", 7, 24, true, 2);
    test_sites(sites, 8, "8mer_BT", 8, 24, true, 2);
    test_sites(sites, 9, "8mer_BT", 9, 24, true, 2);
    test_sites(sites, 10, "8mer_BT", 10, 24, true, 2);

    test_sites(sites, 11, "7mer_BT", 11, 24, true, 3);
    test_sites(sites, 12, "8mer_BT", 12, 24, true, 3);
    test_sites(sites, 13, "8mer_BT", 13, 24, true, 3);
    test_sites(sites, 14, "8mer_BT", 14, 24, true, 3);
    test_sites(sites, 15, "8mer_BT", 15, 24, true, 3);

    test_sites(sites, 16, "7mer_BT", 16, 24, true, 4);
    test_sites(sites, 17, "8mer_BT", 17, 24, true, 4);
    test_sites(sites, 18, "8mer_BT", 18, 24, true, 4);
    test_sites(sites, 19, "8mer_BT", 19, 24, true, 4);

    test_sites(sites, 20, "7mer_BT", 20, 24, true, 5);
    test_sites(sites, 21, "8mer_BT", 21, 24, true, 5);
    test_sites(sites, 22, "8mer_BT", 22, 24, true, 5);
    test_sites(sites, 23, "8mer_BT", 23, 24, true, 5);
    test_sites(sites, 24, "8mer_BT", 24, 24, true, 5);
}

TEST_F(Site05BT1, mir124_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "+";
    mSeedDef[4] = "1:1";
    mSeedDef[5] = "1";
    TSeed seedSeqs;

    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(37u, sites.get_length());

//    test_sites(sites, 0, "MM", 4, 25, false, 0);
//    test_sites(sites, 1, "MM", 5, 25, false, 0);
//
//    test_sites(sites, 2, "MM", 0, 25, false, 0);
//    test_sites(sites, 3, "MM", 1, 25, false, 0);
//    test_sites(sites, 4, "MM", 2, 25, false, 0);
//    test_sites(sites, 5, "MM", 3, 25, false, 0);

    test_sites(sites, 0, "8mer_MM", 9, 24, true, 4);
    test_sites(sites, 1, "8mer_MM", 19, 24, true, 4);

    test_sites(sites, 2, "7mer_MM", 16, 24, true, 4);
    test_sites(sites, 3, "8mer_MM", 17, 24, true, 4);
    test_sites(sites, 4, "8mer_MM", 18, 24, true, 4);

    test_sites(sites, 5, "8mer_MM", 23, 24, true, 5);
    test_sites(sites, 6, "8mer_MM", 24, 24, true, 5);

    test_sites(sites, 7, "7mer_MM", 20, 24, true, 5);
    test_sites(sites, 8, "8mer_MM", 21, 24, true, 5);
    test_sites(sites, 9, "8mer_MM", 22, 24, true, 5);

    test_sites(sites, 10, "8mer_MMGU", 10, 24, true, 4);
    test_sites(sites, 11, "8mer_MMGU", 15, 24, true, 4);

//    test_sites(sites, 18, "BT", 20, 23, false, 0);
//    test_sites(sites, 19, "BT", 21, 23, false, 0);
//    test_sites(sites, 20, "BT", 22, 23, false, 0);
//    test_sites(sites, 21, "BT", 23, 23, false, 0);
//    test_sites(sites, 22, "BT", 24, 23, false, 0);

    test_sites(sites, 12, "7mer_BT", 0, 24, true, 1);
    test_sites(sites, 13, "7mer_BT", 1, 24, true, 1);
    test_sites(sites, 14, "8mer_BT", 2, 24, true, 1);
    test_sites(sites, 15, "8mer_BT", 3, 24, true, 1);
    test_sites(sites, 16, "7mer_BT", 4, 24, true, 1);
    test_sites(sites, 17, "7mer_BT", 5, 24, true, 1);

    test_sites(sites, 18, "7mer_BT", 6, 24, true, 2);
    test_sites(sites, 19, "8mer_BT", 7, 24, true, 2);
    test_sites(sites, 20, "8mer_BT", 8, 24, true, 2);
    test_sites(sites, 21, "8mer_BT", 9, 24, true, 2);
    test_sites(sites, 22, "8mer_BT", 10, 24, true, 2);

    test_sites(sites, 23, "7mer_BT", 11, 24, true, 3);
    test_sites(sites, 24, "8mer_BT", 12, 24, true, 3);
    test_sites(sites, 25, "8mer_BT", 13, 24, true, 3);
    test_sites(sites, 26, "8mer_BT", 14, 24, true, 3);
    test_sites(sites, 27, "8mer_BT", 15, 24, true, 3);

    test_sites(sites, 28, "7mer_BT", 16, 24, true, 4);
    test_sites(sites, 29, "8mer_BT", 17, 24, true, 4);
    test_sites(sites, 30, "8mer_BT", 18, 24, true, 4);
    test_sites(sites, 31, "8mer_BT", 19, 24, true, 4);

    test_sites(sites, 32, "7mer_BT", 20, 24, true, 5);
    test_sites(sites, 33, "8mer_BT", 21, 24, true, 5);
    test_sites(sites, 34, "8mer_BT", 22, 24, true, 5);
    test_sites(sites, 35, "8mer_BT", 23, 24, true, 5);
    test_sites(sites, 36, "8mer_BT", 24, 24, true, 5);
}

}
