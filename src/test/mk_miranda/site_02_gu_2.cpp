#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class Site02GU2 : public TestSiteMR3AS {
protected:
    Site02GU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_miranda/";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "1";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "0";
    }

};

TEST_F(Site02GU2, mir1_gu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    TOp ops;
    TSeed seedSeqs(ops);
    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[1]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(20u, sites.get_length());

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);
    test_sites(sites, 2, "8mer_GUM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_GUM", 3, 24, true, -1);

    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "8mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "8mer_GUT", 7, 24, true, 0);

//    test_sites(sites, 8, "GUT", 8, 24, false, 0);
//    test_sites(sites, 9, "GUT", 9, 24, false, 0);
//    test_sites(sites, 10, "GUT", 10, 24, false, 0);
//    test_sites(sites, 11, "GUT", 11, 24, false, 0);
//    test_sites(sites, 12, "GUT", 12, 24, false, 0);

    test_sites(sites, 8, "7mer_GUM", 13, 24, true, 1);
    test_sites(sites, 9, "7mer_GUM", 14, 24, true, 1);
    test_sites(sites, 10, "8mer_GUM", 15, 24, true, 1);
    test_sites(sites, 11, "8mer_GUM", 16, 24, true, 1);

//    test_sites(sites, 17, "GUM", 17, 24, false, 0);
//    test_sites(sites, 18, "GUM", 18, 24, false, 0);
//    test_sites(sites, 19, "GUM", 19, 24, false, 0);
//    test_sites(sites, 20, "GUM", 20, 24, false, 0);
//    test_sites(sites, 21, "GUM", 21, 24, false, 0);

    test_sites(sites, 12, "7mer_GUT", 22, 24, true, 4);
    test_sites(sites, 13, "7mer_GUT", 23, 24, true, 4);
    test_sites(sites, 14, "8mer_GUT", 24, 24, true, 4);
    test_sites(sites, 15, "8mer_GUT", 25, 24, true, 4);

//    test_sites(sites, 26, "GUT", 26, 24, false, 0);
//    test_sites(sites, 27, "GUT", 27, 24, false, 0);
//    test_sites(sites, 28, "GUT", 28, 24, false, 0);
//    test_sites(sites, 29, "GUT", 29, 24, false, 0);
//    test_sites(sites, 30, "GUT", 30, 24, false, 0);

    test_sites(sites, 16, "7mer_GUT", 31, 24, true, 5);
    test_sites(sites, 17, "7mer_GUT", 32, 24, true, 5);
    test_sites(sites, 18, "8mer_GUT", 33, 24, true, 5);
    test_sites(sites, 19, "8mer_GUT", 34, 24, true, 5);

//    test_sites(sites, 35, "GUT", 35, 24, false, 0);
//    test_sites(sites, 36, "GUT", 36, 24, false, 0);
//    test_sites(sites, 37, "GUT", 37, 24, false, 0);
//    test_sites(sites, 38, "GUT", 38, 24, false, 0);
//    test_sites(sites, 39, "GUT", 39, 24, false, 0);
}

TEST_F(Site02GU2, mir1_gu_plus) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "+";
    TOp ops;
    TSeed seedSeqs(ops);
    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[1]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(36u, sites.get_length());

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);
    test_sites(sites, 2, "8mer_GUM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_GUM", 3, 24, true, -1);

    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "8mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "8mer_GUT", 7, 24, true, 0);

//    test_sites(sites, 8, "GUT", 8, 24, false, 0);
    test_sites(sites, 8, "7mer_GU+", 9, 24, true, 0);
    test_sites(sites, 9, "7mer_GU+", 10, 24, true, 0);
    test_sites(sites, 10, "8mer_GU+", 11, 24, true, 0);
    test_sites(sites, 11, "8mer_GU+", 12, 24, true, 0);

    test_sites(sites, 12, "7mer_GUM", 13, 24, true, 1);
    test_sites(sites, 13, "7mer_GUM", 14, 24, true, 1);
    test_sites(sites, 14, "8mer_GUM", 15, 24, true, 1);
    test_sites(sites, 15, "8mer_GUM", 16, 24, true, 1);

//    test_sites(sites, 17, "GUM", 17, 24, false, 0);
    test_sites(sites, 16, "7mer_GU+", 18, 24, true, 0);
    test_sites(sites, 17, "7mer_GU+", 19, 24, true, 0);
    test_sites(sites, 18, "8mer_GU+", 20, 24, true, 0);
    test_sites(sites, 19, "8mer_GU+", 21, 24, true, 0);

    test_sites(sites, 20, "7mer_GUT", 22, 24, true, 4);
    test_sites(sites, 21, "7mer_GUT", 23, 24, true, 4);
    test_sites(sites, 22, "8mer_GUT", 24, 24, true, 4);
    test_sites(sites, 23, "8mer_GUT", 25, 24, true, 4);

//    test_sites(sites, 26, "GUT", 26, 24, false, 0);
    test_sites(sites, 24, "7mer_GU+", 27, 24, true, 0);
    test_sites(sites, 25, "7mer_GU+", 28, 24, true, 0);
    test_sites(sites, 26, "8mer_GU+", 29, 24, true, 0);
    test_sites(sites, 27, "8mer_GU+", 30, 24, true, 0);

    test_sites(sites, 28, "7mer_GUT", 31, 24, true, 5);
    test_sites(sites, 29, "7mer_GUT", 32, 24, true, 5);
    test_sites(sites, 30, "8mer_GUT", 33, 24, true, 5);
    test_sites(sites, 31, "8mer_GUT", 34, 24, true, 5);

//    test_sites(sites, 35, "GUT", 35, 24, false, 0);
    test_sites(sites, 32, "7mer_GU+", 36, 24, true, 0);
    test_sites(sites, 33, "7mer_GU+", 37, 24, true, 0);
    test_sites(sites, 34, "8mer_GU+", 38, 24, true, 0);
    test_sites(sites, 35, "8mer_GU+", 39, 24, true, 0);
}

TEST_F(Site02GU2, mir1_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "+";
    mSeedDef[4] = "1:1";
    mSeedDef[5] = "1";
    TOp ops;
    TSeed seedSeqs(ops);
    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[1]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(83u, sites.get_length());

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);
    test_sites(sites, 2, "8mer_GUM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_GUM", 3, 24, true, -1);

    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "8mer_GUT", 6, 24, true, 0);
    test_sites(sites, 7, "8mer_GUT", 7, 24, true, 0);

    test_sites(sites, 8, "7mer_MMGU", 8, 24, true, -1);
    test_sites(sites, 9, "7mer_GU+", 9, 24, true, 0);
    test_sites(sites, 10, "7mer_GU+", 10, 24, true, 0);
    test_sites(sites, 11, "8mer_GU+", 11, 24, true, 0);
    test_sites(sites, 12, "8mer_GU+", 12, 24, true, 0);

    test_sites(sites, 13, "7mer_GUM", 13, 24, true, 1);
    test_sites(sites, 14, "7mer_GUM", 14, 24, true, 1);
    test_sites(sites, 15, "8mer_GUM", 15, 24, true, 1);
    test_sites(sites, 16, "8mer_GUM", 16, 24, true, 1);

    test_sites(sites, 17, "7mer_MMGU", 17, 24, true, -1);
    test_sites(sites, 18, "7mer_GU+", 18, 24, true, 0);
    test_sites(sites, 19, "7mer_GU+", 19, 24, true, 0);
    test_sites(sites, 20, "8mer_GU+", 20, 24, true, 0);
    test_sites(sites, 21, "8mer_GU+", 21, 24, true, 0);

    test_sites(sites, 22, "7mer_GUT", 22, 24, true, 4);
    test_sites(sites, 23, "7mer_GUT", 23, 24, true, 4);
    test_sites(sites, 24, "8mer_GUT", 24, 24, true, 4);
    test_sites(sites, 25, "8mer_GUT", 25, 24, true, 4);

    test_sites(sites, 26, "7mer_MMGU", 26, 24, true, -1);
    test_sites(sites, 27, "7mer_GU+", 27, 24, true, 0);
    test_sites(sites, 28, "7mer_GU+", 28, 24, true, 0);
    test_sites(sites, 29, "8mer_GU+", 29, 24, true, 0);
    test_sites(sites, 30, "8mer_GU+", 30, 24, true, 0);

    test_sites(sites, 31, "7mer_GUT", 31, 24, true, 5);
    test_sites(sites, 32, "7mer_GUT", 32, 24, true, 5);
    test_sites(sites, 33, "8mer_GUT", 33, 24, true, 5);
    test_sites(sites, 34, "8mer_GUT", 34, 24, true, 5);

    test_sites(sites, 35, "7mer_MMGU", 35, 24, true, -1);
    test_sites(sites, 36, "7mer_GU+", 36, 24, true, 0);
    test_sites(sites, 37, "7mer_GU+", 37, 24, true, 0);
    test_sites(sites, 38, "8mer_GU+", 38, 24, true, 0);
    test_sites(sites, 39, "8mer_GU+", 39, 24, true, 0);

    test_sites(sites, 40, "7mer_MMGU", 0, 12, true, 3);
    test_sites(sites, 41, "7mer_MMGU", 1, 12, true, 3);
    test_sites(sites, 42, "7mer_MMGU", 2, 12, true, 3);
    test_sites(sites, 43, "7mer_MMGU", 3, 12, true, 3);

    test_sites(sites, 44, "7mer_MMGU", 4, 12, true, 3);
    test_sites(sites, 45, "7mer_MMGU", 5, 12, true, 3);
    test_sites(sites, 46, "7mer_MMGU", 6, 12, true, 3);
    test_sites(sites, 47, "7mer_MMGU", 7, 12, true, 3);

    test_sites(sites, 48, "7mer_MMGU", 8, 12, true, 3);
    test_sites(sites, 49, "7mer_MMGU", 9, 12, true, 3);
    test_sites(sites, 50, "7mer_MMGU", 10, 12, true, 3);
    test_sites(sites, 51, "7mer_MMGU", 11, 12, true, 3);
    test_sites(sites, 52, "7mer_MMGU", 12, 12, true, 3);

    test_sites(sites, 53, "7mer_MMGU", 13, 12, true, 3);
    test_sites(sites, 54, "7mer_MMGU", 14, 12, true, 3);
    test_sites(sites, 55, "7mer_MMGU", 15, 12, true, 3);
    test_sites(sites, 56, "7mer_MMGU", 16, 12, true, 3);

    test_sites(sites, 57, "7mer_MMGU", 17, 12, true, 3);
    test_sites(sites, 58, "7mer_MMGU", 18, 12, true, 3);
    test_sites(sites, 59, "7mer_MMGU", 19, 12, true, 3);
    test_sites(sites, 60, "7mer_MMGU", 20, 12, true, 3);
    test_sites(sites, 61, "7mer_MMGU", 21, 12, true, 3);

    test_sites(sites, 62, "7mer_MMGU", 22, 12, true, 3);
    test_sites(sites, 63, "7mer_MMGU", 23, 12, true, 3);
    test_sites(sites, 64, "7mer_MMGU", 24, 12, true, 3);
    test_sites(sites, 65, "7mer_MMGU", 25, 12, true, 3);

    test_sites(sites, 66, "7mer_MMGU", 26, 12, true, 3);
    test_sites(sites, 67, "7mer_MMGU", 27, 12, true, 3);
    test_sites(sites, 68, "7mer_MMGU", 28, 12, true, 3);
    test_sites(sites, 69, "7mer_MMGU", 29, 12, true, 3);
    test_sites(sites, 70, "7mer_MMGU", 30, 12, true, 3);

    test_sites(sites, 71, "7mer_MMGU", 31, 12, true, 3);
    test_sites(sites, 72, "7mer_MMGU", 32, 12, true, 3);
    test_sites(sites, 73, "7mer_MMGU", 33, 12, true, 3);
    test_sites(sites, 74, "7mer_MMGU", 34, 12, true, 3);

    test_sites(sites, 75, "7mer_MMGU", 35, 12, true, 3);
    test_sites(sites, 76, "7mer_MMGU", 36, 12, true, 3);
    test_sites(sites, 77, "7mer_MMGU", 37, 12, true, 3);
    test_sites(sites, 78, "7mer_MMGU", 38, 12, true, 3);
    test_sites(sites, 79, "7mer_MMGU", 39, 12, true, 3);

//    test_sites(sites, 80, "BT", 31, 23, false, 0);
//    test_sites(sites, 81, "BT", 32, 23, false, 0);
//    test_sites(sites, 82, "BT", 33, 23, false, 0);
//    test_sites(sites, 83, "BT", 34, 23, false, 0);
//    test_sites(sites, 84, "BT", 35, 23, false, 0);

    test_sites(sites, 80, "7mer_BT", 0, 23, true, 0);
    test_sites(sites, 81, "7mer_BT", 1, 23, true, 0);
//    test_sites(sites, 87, "BT", 2, 23, false, 0);
//    test_sites(sites, 88, "BT", 3, 23, false, 0);

//    test_sites(sites, 89, "BT", 36, 23, false, 0);
//    test_sites(sites, 90, "BT", 37, 23, false, 0);
//    test_sites(sites, 91, "BT", 38, 23, false, 0);
//    test_sites(sites, 92, "BT", 39, 23, false, 0);

    test_sites(sites, 82, "7mer_BT", 8, 23, true, 1);
}
}
