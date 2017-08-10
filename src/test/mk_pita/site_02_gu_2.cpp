#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site02GU2 : public TestSitePITA {
protected:
    Site02GU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_pita/";

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

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
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

    TOp ops;
    TSeed seedSeqs(ops);
    mSeedDef[3] = "+";
    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[1]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(28u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
    test_sites(sites, 2, "8mer_GUM", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_GUM", 3, 24, true, -1);
    test_sites(sites, 4, "7mer_GU+", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GU+", 5, 24, true, 0);
    test_sites(sites, 6, "8mer_GU+", 6, 24, true, 0);
    test_sites(sites, 7, "8mer_GU+", 7, 24, true, 0);
//    test_sites(sites, 8, "GU+", 8, 24, false, 0);
//    test_sites(sites, 9, "GU+", 9, 24, false, 0);
//    test_sites(sites, 10, "GU+", 10, 24, false, 0);
    test_sites(sites, 8, "8mer_GU+", 11, 24, true, 0);
    test_sites(sites, 9, "8mer_GU+", 12, 24, true, 0);
    test_sites(sites, 10, "7mer_GU+", 13, 24, true, 0);
    test_sites(sites, 11, "7mer_GU+", 14, 24, true, 0);
    test_sites(sites, 12, "8mer_GU+", 15, 24, true, 0);
    test_sites(sites, 13, "8mer_GU+", 16, 24, true, 0);
//    test_sites(sites, 17, "GU+", 17, 24, false, 0);
//    test_sites(sites, 18, "GU+", 18, 24, false, 0);
//    test_sites(sites, 19, "GU+", 19, 24, false, 0);
    test_sites(sites, 14, "8mer_GU+", 20, 24, true, 0);
    test_sites(sites, 15, "8mer_GU+", 21, 24, true, 0);
    test_sites(sites, 16, "7mer_GU+", 22, 24, true, 0);
    test_sites(sites, 17, "7mer_GU+", 23, 24, true, 0);
    test_sites(sites, 18, "8mer_GU+", 24, 24, true, 0);
    test_sites(sites, 19, "8mer_GU+", 25, 24, true, 0);
//    test_sites(sites, 26, "GU+", 26, 24, false, 0);
//    test_sites(sites, 27, "GU+", 27, 24, false, 0);
//    test_sites(sites, 28, "GU+", 28, 24, false, 0);
    test_sites(sites, 20, "8mer_GU+", 29, 24, true, 0);
    test_sites(sites, 21, "8mer_GU+", 30, 24, true, 0);
    test_sites(sites, 22, "7mer_GU+", 31, 24, true, 0);
    test_sites(sites, 23, "7mer_GU+", 32, 24, true, 0);
    test_sites(sites, 24, "8mer_GU+", 33, 24, true, 0);
    test_sites(sites, 25, "8mer_GU+", 34, 24, true, 0);
//    test_sites(sites, 35, "GU+", 35, 24, false, 0);
//    test_sites(sites, 36, "GU+", 36, 24, false, 0);
//    test_sites(sites, 37, "GU+", 37, 24, false, 0);
    test_sites(sites, 26, "8mer_GU+", 38, 24, true, 0);
    test_sites(sites, 27, "8mer_GU+", 39, 24, true, 0);

}

TEST_F(Site02GU2, mir1_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    TOp ops;
    TSeed seedSeqs(ops);
    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[1]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(20u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);
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
//    test_sites(sites, 40, "", 0, 12, false, 0);
//    test_sites(sites, 41, "", 1, 12, false, 0);
//    test_sites(sites, 42, "", 2, 12, false, 0);
//    test_sites(sites, 43, "", 3, 12, false, 0);
//    test_sites(sites, 44, "", 4, 12, false, 0);
//    test_sites(sites, 45, "", 5, 12, false, 0);
//    test_sites(sites, 46, "", 6, 12, false, 0);
//    test_sites(sites, 47, "", 7, 12, false, 0);
//    test_sites(sites, 48, "", 8, 12, false, 0);
//    test_sites(sites, 49, "", 9, 12, false, 0);
//    test_sites(sites, 50, "", 10, 12, false, 0);
//    test_sites(sites, 51, "", 11, 12, false, 0);
//    test_sites(sites, 52, "", 12, 12, false, 0);
//    test_sites(sites, 53, "", 13, 12, false, 0);
//    test_sites(sites, 54, "", 14, 12, false, 0);
//    test_sites(sites, 55, "", 15, 12, false, 0);
//    test_sites(sites, 56, "", 16, 12, false, 0);
//    test_sites(sites, 57, "", 17, 12, false, 0);
//    test_sites(sites, 58, "", 18, 12, false, 0);
//    test_sites(sites, 59, "", 19, 12, false, 0);
//    test_sites(sites, 60, "", 20, 12, false, 0);
//    test_sites(sites, 61, "", 21, 12, false, 0);
//    test_sites(sites, 62, "", 22, 12, false, 0);
//    test_sites(sites, 63, "", 23, 12, false, 0);
//    test_sites(sites, 64, "", 24, 12, false, 0);
//    test_sites(sites, 65, "", 25, 12, false, 0);
//    test_sites(sites, 66, "", 26, 12, false, 0);
//    test_sites(sites, 67, "", 27, 12, false, 0);
//    test_sites(sites, 68, "", 28, 12, false, 0);
//    test_sites(sites, 69, "", 29, 12, false, 0);
//    test_sites(sites, 70, "", 30, 12, false, 0);
//    test_sites(sites, 71, "", 31, 12, false, 0);
//    test_sites(sites, 72, "", 32, 12, false, 0);
//    test_sites(sites, 73, "", 33, 12, false, 0);
//    test_sites(sites, 74, "", 34, 12, false, 0);
//    test_sites(sites, 75, "", 35, 12, false, 0);
//    test_sites(sites, 76, "", 36, 12, false, 0);
//    test_sites(sites, 77, "", 37, 12, false, 0);
//    test_sites(sites, 78, "", 38, 12, false, 0);
//    test_sites(sites, 79, "", 39, 12, false, 0);

}

}
