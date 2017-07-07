#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site04MMGU2 : public TestSitePITA {
protected:
    Site04MMGU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_2.fasta";
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
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "0";
    }

    typedef ptddg::PITACore<mikan::TRNATYPE>::TIndexQGram TIdx;
    typedef ptddg::PITACore<mikan::TRNATYPE>::TFinder TFin;
    typedef ptddg::PITASeedSites<mikan::TRNATYPE> TSit;

};

TEST_F(Site04MMGU2, mir1_mm7gu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(94u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);

    test_sites(sites, 2, "8mer_MMGU", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MMGU", 3, 24, true, -1);
    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "GUT", 6, 24, false, 0);
    test_sites(sites, 7, "GUT", 7, 24, false, 0);

    test_sites(sites, 8, "8mer_MMGU", 8, 24, true, -1);
    test_sites(sites, 9, "8mer_MMGU", 9, 24, true, -1);
    test_sites(sites, 10, "7mer_GUM", 10, 24, true, 1);
    test_sites(sites, 11, "7mer_GUM", 11, 24, true, 1);
    test_sites(sites, 12, "GUM", 12, 24, false, 0);
    test_sites(sites, 13, "GUM", 13, 24, false, 0);

    test_sites(sites, 14, "8mer_MMGU", 14, 24, true, -1);
    test_sites(sites, 15, "8mer_MMGU", 15, 24, true, -1);
    test_sites(sites, 16, "7mer_GUT", 16, 24, true, 4);
    test_sites(sites, 17, "7mer_GUT", 17, 24, true, 4);
    test_sites(sites, 18, "GUT", 18, 24, false, 0);
    test_sites(sites, 19, "GUT", 19, 24, false, 0);

    test_sites(sites, 20, "8mer_MMGU", 20, 24, true, -1);
    test_sites(sites, 21, "8mer_MMGU", 21, 24, true, -1);
    test_sites(sites, 22, "7mer_GUT", 22, 24, true, 5);
    test_sites(sites, 23, "7mer_GUT", 23, 24, true, 5);
    test_sites(sites, 24, "GUT", 24, 24, false, 0);
    test_sites(sites, 25, "GUT", 25, 24, false, 0);

    test_sites(sites, 26, "8mer_MMGU", 26, 24, true, 0);
    test_sites(sites, 27, "8mer_MMGU", 27, 24, true, 0);
    test_sites(sites, 28, "8mer_MMGU", 28, 24, true, 1);
    test_sites(sites, 29, "8mer_MMGU", 29, 24, true, 1);
    test_sites(sites, 30, "8mer_MMGU", 30, 24, true, 2);
    test_sites(sites, 31, "8mer_MMGU", 31, 24, true, 2);
    test_sites(sites, 32, "8mer_MMGU", 32, 24, true, 3);
    test_sites(sites, 33, "8mer_MMGU", 33, 24, true, 3);
    test_sites(sites, 34, "8mer_MMGU", 34, 24, true, 4);
    test_sites(sites, 35, "8mer_MMGU", 35, 24, true, 4);
    test_sites(sites, 36, "8mer_MMGU", 36, 24, true, 5);
    test_sites(sites, 37, "8mer_MMGU", 37, 24, true, 5);

    test_sites(sites, 38, "8mer_MMGU", 38, 24, true, 1);
    test_sites(sites, 39, "8mer_MMGU", 39, 24, true, 1);
    test_sites(sites, 40, "7mer_MMGU", 40, 24, true, 1);
    test_sites(sites, 41, "7mer_MMGU", 41, 24, true, 1);
    test_sites(sites, 42, "8mer_MMGU", 42, 24, true, 2);
    test_sites(sites, 43, "8mer_MMGU", 43, 24, true, 2);
    test_sites(sites, 44, "7mer_MMGU", 44, 24, true, 2);
    test_sites(sites, 45, "7mer_MMGU", 45, 24, true, 2);
    test_sites(sites, 46, "8mer_MMGU", 46, 24, true, 3);
    test_sites(sites, 47, "8mer_MMGU", 47, 24, true, 3);
    test_sites(sites, 48, "7mer_MMGU", 48, 24, true, 3);
    test_sites(sites, 49, "7mer_MMGU", 49, 24, true, 3);
    test_sites(sites, 50, "8mer_MMGU", 50, 24, true, 4);
    test_sites(sites, 51, "8mer_MMGU", 51, 24, true, 4);
    test_sites(sites, 52, "7mer_MMGU", 52, 24, true, 4);
    test_sites(sites, 53, "7mer_MMGU", 53, 24, true, 4);
    test_sites(sites, 54, "8mer_MMGU", 54, 24, true, 5);
    test_sites(sites, 55, "8mer_MMGU", 55, 24, true, 5);
    test_sites(sites, 56, "7mer_MMGU", 56, 24, true, 5);
    test_sites(sites, 57, "7mer_MMGU", 57, 24, true, 5);

    test_sites(sites, 58, "8mer_MMGU", 58, 24, true, 0);
    test_sites(sites, 59, "8mer_MMGU", 59, 24, true, 0);
    test_sites(sites, 60, "7mer_MMGU", 60, 24, true, 0);
    test_sites(sites, 61, "7mer_MMGU", 61, 24, true, 0);
    test_sites(sites, 62, "8mer_MMGU", 62, 24, true, 2);
    test_sites(sites, 63, "8mer_MMGU", 63, 24, true, 2);
    test_sites(sites, 64, "7mer_MMGU", 64, 24, true, 2);
    test_sites(sites, 65, "7mer_MMGU", 65, 24, true, 2);
    test_sites(sites, 66, "8mer_MMGU", 66, 24, true, 3);
    test_sites(sites, 67, "8mer_MMGU", 67, 24, true, 3);
    test_sites(sites, 68, "7mer_MMGU", 68, 24, true, 3);
    test_sites(sites, 69, "7mer_MMGU", 69, 24, true, 3);
    test_sites(sites, 70, "8mer_MMGU", 70, 24, true, 4);
    test_sites(sites, 71, "8mer_MMGU", 71, 24, true, 4);
    test_sites(sites, 72, "7mer_MMGU", 72, 24, true, 4);
    test_sites(sites, 73, "7mer_MMGU", 73, 24, true, 4);
    test_sites(sites, 74, "8mer_MMGU", 74, 24, true, 5);
    test_sites(sites, 75, "8mer_MMGU", 75, 24, true, 5);
    test_sites(sites, 76, "7mer_MMGU", 76, 24, true, 5);
    test_sites(sites, 77, "7mer_MMGU", 77, 24, true, 5);

    test_sites(sites, 78, "8mer_MMGU", 78, 24, true, 0);
    test_sites(sites, 79, "8mer_MMGU", 79, 24, true, 0);
    test_sites(sites, 80, "7mer_MMGU", 80, 24, true, 0);
    test_sites(sites, 81, "7mer_MMGU", 81, 24, true, 0);
    test_sites(sites, 82, "8mer_MMGU", 82, 24, true, 2);
    test_sites(sites, 83, "8mer_MMGU", 83, 24, true, 2);
    test_sites(sites, 84, "7mer_MMGU", 84, 24, true, 2);
    test_sites(sites, 85, "7mer_MMGU", 85, 24, true, 2);
    test_sites(sites, 86, "8mer_MMGU", 86, 24, true, 3);
    test_sites(sites, 87, "8mer_MMGU", 87, 24, true, 3);
    test_sites(sites, 88, "7mer_MMGU", 88, 24, true, 3);
    test_sites(sites, 89, "7mer_MMGU", 89, 24, true, 3);
    test_sites(sites, 90, "8mer_MMGU", 90, 24, true, 5);
    test_sites(sites, 91, "8mer_MMGU", 91, 24, true, 5);
    test_sites(sites, 92, "7mer_MMGU", 92, 24, true, 5);
    test_sites(sites, 93, "7mer_MMGU", 93, 24, true, 5);
}

TEST_F(Site04MMGU2, mir1_mm8gu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[4] = "0:1";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(94u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);

    test_sites(sites, 2, "8mer_MMGU", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MMGU", 3, 24, true, -1);
    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "GUT", 6, 24, false, 0);
    test_sites(sites, 7, "GUT", 7, 24, false, 0);

    test_sites(sites, 8, "8mer_MMGU", 8, 24, true, -1);
    test_sites(sites, 9, "8mer_MMGU", 9, 24, true, -1);
    test_sites(sites, 10, "7mer_GUM", 10, 24, true, 1);
    test_sites(sites, 11, "7mer_GUM", 11, 24, true, 1);
    test_sites(sites, 12, "GUM", 12, 24, false, 0);
    test_sites(sites, 13, "GUM", 13, 24, false, 0);

    test_sites(sites, 14, "8mer_MMGU", 14, 24, true, -1);
    test_sites(sites, 15, "8mer_MMGU", 15, 24, true, -1);
    test_sites(sites, 16, "7mer_GUT", 16, 24, true, 4);
    test_sites(sites, 17, "7mer_GUT", 17, 24, true, 4);
    test_sites(sites, 18, "GUT", 18, 24, false, 0);
    test_sites(sites, 19, "GUT", 19, 24, false, 0);

    test_sites(sites, 20, "8mer_MMGU", 20, 24, true, -1);
    test_sites(sites, 21, "8mer_MMGU", 21, 24, true, -1);
    test_sites(sites, 22, "7mer_GUT", 22, 24, true, 5);
    test_sites(sites, 23, "7mer_GUT", 23, 24, true, 5);
    test_sites(sites, 24, "GUT", 24, 24, false, 0);
    test_sites(sites, 25, "GUT", 25, 24, false, 0);

    test_sites(sites, 26, "8mer_MMGU", 26, 24, true, 0);
    test_sites(sites, 27, "8mer_MMGU", 27, 24, true, 0);
    test_sites(sites, 28, "8mer_MMGU", 28, 24, true, 1);
    test_sites(sites, 29, "8mer_MMGU", 29, 24, true, 1);
    test_sites(sites, 30, "8mer_MMGU", 30, 24, true, 2);
    test_sites(sites, 31, "8mer_MMGU", 31, 24, true, 2);
    test_sites(sites, 32, "8mer_MMGU", 32, 24, true, 3);
    test_sites(sites, 33, "8mer_MMGU", 33, 24, true, 3);
    test_sites(sites, 34, "8mer_MMGU", 34, 24, true, 4);
    test_sites(sites, 35, "8mer_MMGU", 35, 24, true, 4);
    test_sites(sites, 36, "8mer_MMGU", 36, 24, true, 5);
    test_sites(sites, 37, "8mer_MMGU", 37, 24, true, 5);

    test_sites(sites, 38, "8mer_MMGU", 38, 24, true, 1);
    test_sites(sites, 39, "8mer_MMGU", 39, 24, true, 1);
    test_sites(sites, 40, "MMGU", 40, 24, false, 0);
    test_sites(sites, 41, "MMGU", 41, 24, false, 0);
    test_sites(sites, 42, "8mer_MMGU", 42, 24, true, 2);
    test_sites(sites, 43, "8mer_MMGU", 43, 24, true, 2);
    test_sites(sites, 44, "MMGU", 44, 24, false, 0);
    test_sites(sites, 45, "MMGU", 45, 24, false, 0);
    test_sites(sites, 46, "8mer_MMGU", 46, 24, true, 3);
    test_sites(sites, 47, "8mer_MMGU", 47, 24, true, 3);
    test_sites(sites, 48, "MMGU", 48, 24, false, 0);
    test_sites(sites, 49, "MMGU", 49, 24, false, 0);
    test_sites(sites, 50, "8mer_MMGU", 50, 24, true, 4);
    test_sites(sites, 51, "8mer_MMGU", 51, 24, true, 4);
    test_sites(sites, 52, "MMGU", 52, 24, false, 0);
    test_sites(sites, 53, "MMGU", 53, 24, false, 0);
    test_sites(sites, 54, "8mer_MMGU", 54, 24, true, 5);
    test_sites(sites, 55, "8mer_MMGU", 55, 24, true, 5);
    test_sites(sites, 56, "MMGU", 56, 24, false, 0);
    test_sites(sites, 57, "MMGU", 57, 24, false, 0);

    test_sites(sites, 58, "8mer_MMGU", 58, 24, true, 0);
    test_sites(sites, 59, "8mer_MMGU", 59, 24, true, 0);
    test_sites(sites, 60, "MMGU", 60, 24, false, 0);
    test_sites(sites, 61, "MMGU", 61, 24, false, 0);
    test_sites(sites, 62, "8mer_MMGU", 62, 24, true, 2);
    test_sites(sites, 63, "8mer_MMGU", 63, 24, true, 2);
    test_sites(sites, 64, "MMGU", 64, 24, false, 0);
    test_sites(sites, 65, "MMGU", 65, 24, false, 0);
    test_sites(sites, 66, "8mer_MMGU", 66, 24, true, 3);
    test_sites(sites, 67, "8mer_MMGU", 67, 24, true, 3);
    test_sites(sites, 68, "MMGU", 68, 24, false, 0);
    test_sites(sites, 69, "MMGU", 69, 24, false, 0);
    test_sites(sites, 70, "8mer_MMGU", 70, 24, true, 4);
    test_sites(sites, 71, "8mer_MMGU", 71, 24, true, 4);
    test_sites(sites, 72, "MMGU", 72, 24, false, 0);
    test_sites(sites, 73, "MMGU", 73, 24, false, 0);
    test_sites(sites, 74, "8mer_MMGU", 74, 24, true, 5);
    test_sites(sites, 75, "8mer_MMGU", 75, 24, true, 5);
    test_sites(sites, 76, "MMGU", 76, 24, false, 0);
    test_sites(sites, 77, "MMGU", 77, 24, false, 0);

    test_sites(sites, 78, "8mer_MMGU", 78, 24, true, 0);
    test_sites(sites, 79, "8mer_MMGU", 79, 24, true, 0);
    test_sites(sites, 80, "MMGU", 80, 24, false, 0);
    test_sites(sites, 81, "MMGU", 81, 24, false, 0);
    test_sites(sites, 82, "8mer_MMGU", 82, 24, true, 2);
    test_sites(sites, 83, "8mer_MMGU", 83, 24, true, 2);
    test_sites(sites, 84, "MMGU", 84, 24, false, 0);
    test_sites(sites, 85, "MMGU", 85, 24, false, 0);
    test_sites(sites, 86, "8mer_MMGU", 86, 24, true, 3);
    test_sites(sites, 87, "8mer_MMGU", 87, 24, true, 3);
    test_sites(sites, 88, "MMGU", 88, 24, false, 0);
    test_sites(sites, 89, "MMGU", 89, 24, false, 0);
    test_sites(sites, 90, "8mer_MMGU", 90, 24, true, 5);
    test_sites(sites, 91, "8mer_MMGU", 91, 24, true, 5);
    test_sites(sites, 92, "MMGU", 92, 24, false, 0);
    test_sites(sites, 93, "MMGU", 93, 24, false, 0);
}

TEST_F(Site04MMGU2, mir1_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(94u, sites.get_length());

    test_sites(sites, 0, "6mer", 0, 24, true, 0);
    test_sites(sites, 1, "6mer", 1, 24, true, 0);

    test_sites(sites, 2, "8mer_MMGU", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MMGU", 3, 24, true, -1);
    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "GUT", 6, 24, false, 0);
    test_sites(sites, 7, "GUT", 7, 24, false, 0);

    test_sites(sites, 8, "8mer_MMGU", 8, 24, true, -1);
    test_sites(sites, 9, "8mer_MMGU", 9, 24, true, -1);
    test_sites(sites, 10, "7mer_GUM", 10, 24, true, 1);
    test_sites(sites, 11, "7mer_GUM", 11, 24, true, 1);
    test_sites(sites, 12, "GUM", 12, 24, false, 0);
    test_sites(sites, 13, "GUM", 13, 24, false, 0);

    test_sites(sites, 14, "8mer_MMGU", 14, 24, true, -1);
    test_sites(sites, 15, "8mer_MMGU", 15, 24, true, -1);
    test_sites(sites, 16, "7mer_GUT", 16, 24, true, 4);
    test_sites(sites, 17, "7mer_GUT", 17, 24, true, 4);
    test_sites(sites, 18, "GUT", 18, 24, false, 0);
    test_sites(sites, 19, "GUT", 19, 24, false, 0);

    test_sites(sites, 20, "8mer_MMGU", 20, 24, true, -1);
    test_sites(sites, 21, "8mer_MMGU", 21, 24, true, -1);
    test_sites(sites, 22, "7mer_GUT", 22, 24, true, 5);
    test_sites(sites, 23, "7mer_GUT", 23, 24, true, 5);
    test_sites(sites, 24, "GUT", 24, 24, false, 0);
    test_sites(sites, 25, "GUT", 25, 24, false, 0);

    test_sites(sites, 26, "8mer_MMGU", 26, 24, true, 0);
    test_sites(sites, 27, "8mer_MMGU", 27, 24, true, 0);
    test_sites(sites, 28, "8mer_MMGU", 28, 24, true, 1);
    test_sites(sites, 29, "8mer_MMGU", 29, 24, true, 1);
    test_sites(sites, 30, "8mer_MMGU", 30, 24, true, 2);
    test_sites(sites, 31, "8mer_MMGU", 31, 24, true, 2);
    test_sites(sites, 32, "8mer_MMGU", 32, 24, true, 3);
    test_sites(sites, 33, "8mer_MMGU", 33, 24, true, 3);
    test_sites(sites, 34, "8mer_MMGU", 34, 24, true, 4);
    test_sites(sites, 35, "8mer_MMGU", 35, 24, true, 4);
    test_sites(sites, 36, "8mer_MMGU", 36, 24, true, 5);
    test_sites(sites, 37, "8mer_MMGU", 37, 24, true, 5);

    test_sites(sites, 38, "8mer_MMGU", 38, 24, true, 1);
    test_sites(sites, 39, "8mer_MMGU", 39, 24, true, 1);
    test_sites(sites, 40, "MMGU", 40, 24, false, 0);
    test_sites(sites, 41, "MMGU", 41, 24, false, 0);
    test_sites(sites, 42, "8mer_MMGU", 42, 24, true, 2);
    test_sites(sites, 43, "8mer_MMGU", 43, 24, true, 2);
    test_sites(sites, 44, "MMGU", 44, 24, false, 0);
    test_sites(sites, 45, "MMGU", 45, 24, false, 0);
    test_sites(sites, 46, "8mer_MMGU", 46, 24, true, 3);
    test_sites(sites, 47, "8mer_MMGU", 47, 24, true, 3);
    test_sites(sites, 48, "MMGU", 48, 24, false, 0);
    test_sites(sites, 49, "MMGU", 49, 24, false, 0);
    test_sites(sites, 50, "8mer_MMGU", 50, 24, true, 4);
    test_sites(sites, 51, "8mer_MMGU", 51, 24, true, 4);
    test_sites(sites, 52, "MMGU", 52, 24, false, 0);
    test_sites(sites, 53, "MMGU", 53, 24, false, 0);
    test_sites(sites, 54, "8mer_MMGU", 54, 24, true, 5);
    test_sites(sites, 55, "8mer_MMGU", 55, 24, true, 5);
    test_sites(sites, 56, "MMGU", 56, 24, false, 0);
    test_sites(sites, 57, "MMGU", 57, 24, false, 0);

    test_sites(sites, 58, "8mer_MMGU", 58, 24, true, 0);
    test_sites(sites, 59, "8mer_MMGU", 59, 24, true, 0);
    test_sites(sites, 60, "MMGU", 60, 24, false, 0);
    test_sites(sites, 61, "MMGU", 61, 24, false, 0);
    test_sites(sites, 62, "8mer_MMGU", 62, 24, true, 2);
    test_sites(sites, 63, "8mer_MMGU", 63, 24, true, 2);
    test_sites(sites, 64, "MMGU", 64, 24, false, 0);
    test_sites(sites, 65, "MMGU", 65, 24, false, 0);
    test_sites(sites, 66, "8mer_MMGU", 66, 24, true, 3);
    test_sites(sites, 67, "8mer_MMGU", 67, 24, true, 3);
    test_sites(sites, 68, "MMGU", 68, 24, false, 0);
    test_sites(sites, 69, "MMGU", 69, 24, false, 0);
    test_sites(sites, 70, "8mer_MMGU", 70, 24, true, 4);
    test_sites(sites, 71, "8mer_MMGU", 71, 24, true, 4);
    test_sites(sites, 72, "MMGU", 72, 24, false, 0);
    test_sites(sites, 73, "MMGU", 73, 24, false, 0);
    test_sites(sites, 74, "8mer_MMGU", 74, 24, true, 5);
    test_sites(sites, 75, "8mer_MMGU", 75, 24, true, 5);
    test_sites(sites, 76, "MMGU", 76, 24, false, 0);
    test_sites(sites, 77, "MMGU", 77, 24, false, 0);

    test_sites(sites, 78, "8mer_MMGU", 78, 24, true, 0);
    test_sites(sites, 79, "8mer_MMGU", 79, 24, true, 0);
    test_sites(sites, 80, "MMGU", 80, 24, false, 0);
    test_sites(sites, 81, "MMGU", 81, 24, false, 0);
    test_sites(sites, 82, "8mer_MMGU", 82, 24, true, 2);
    test_sites(sites, 83, "8mer_MMGU", 83, 24, true, 2);
    test_sites(sites, 84, "MMGU", 84, 24, false, 0);
    test_sites(sites, 85, "MMGU", 85, 24, false, 0);
    test_sites(sites, 86, "8mer_MMGU", 86, 24, true, 3);
    test_sites(sites, 87, "8mer_MMGU", 87, 24, true, 3);
    test_sites(sites, 88, "MMGU", 88, 24, false, 0);
    test_sites(sites, 89, "MMGU", 89, 24, false, 0);
    test_sites(sites, 90, "8mer_MMGU", 90, 24, true, 5);
    test_sites(sites, 91, "8mer_MMGU", 91, 24, true, 5);
    test_sites(sites, 92, "MMGU", 92, 24, false, 0);
    test_sites(sites, 93, "MMGU", 93, 24, false, 0);
}
}
