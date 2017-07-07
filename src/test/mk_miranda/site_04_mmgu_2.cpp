#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class Site04MMGU2 : public TestSiteMR3AS {
protected:
    Site04MMGU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_2.fasta";
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
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "0";
    }

    typedef mr3as::MR3Core<mikan::TRNATYPE>::TIndexQGram TIdx;
    typedef mr3as::MR3Core<mikan::TRNATYPE>::TFinder TFin;
    typedef mr3as::MR3SeedSites<mikan::TRNATYPE> TSit;

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

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);

    test_sites(sites, 2, "8mer_MMGU", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MMGU", 3, 24, true, -1);
    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "7mer_MMGU", 6, 24, true, -1);
    test_sites(sites, 7, "7mer_MMGU", 7, 24, true, -1);

    test_sites(sites, 8, "8mer_MMGU", 8, 24, true, -1);
    test_sites(sites, 9, "8mer_MMGU", 9, 24, true, -1);
    test_sites(sites, 10, "7mer_GUM", 10, 24, true, 1);
    test_sites(sites, 11, "7mer_GUM", 11, 24, true, 1);
    test_sites(sites, 12, "7mer_MMGU", 12, 24, true, -1);
    test_sites(sites, 13, "7mer_MMGU", 13, 24, true, -1);

    test_sites(sites, 14, "8mer_MMGU", 14, 24, true, -1);
    test_sites(sites, 15, "8mer_MMGU", 15, 24, true, -1);
    test_sites(sites, 16, "7mer_GUT", 16, 24, true, 4);
    test_sites(sites, 17, "7mer_GUT", 17, 24, true, 4);
    test_sites(sites, 18, "7mer_MMGU", 18, 24, true, -1);
    test_sites(sites, 19, "7mer_MMGU", 19, 24, true, -1);

    test_sites(sites, 20, "8mer_MMGU", 20, 24, true, -1);
    test_sites(sites, 21, "8mer_MMGU", 21, 24, true, -1);
    test_sites(sites, 22, "7mer_GUT", 22, 24, true, 5);
    test_sites(sites, 23, "7mer_GUT", 23, 24, true, 5);
    test_sites(sites, 24, "7mer_MMGU", 24, 24, true, -1);
    test_sites(sites, 25, "7mer_MMGU", 25, 24, true, -1);

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

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);

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

    mSeedDef[3] = "+";
    mSeedDef[4] = "1:1";
    mSeedDef[5] = "1";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(212u, sites.get_length());

    test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
    test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);
    test_sites(sites, 2, "8mer_MMGU", 2, 24, true, -1);
    test_sites(sites, 3, "8mer_MMGU", 3, 24, true, -1);
    test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
    test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
    test_sites(sites, 6, "7mer_MMGU", 6, 24, true, -1);
    test_sites(sites, 7, "7mer_MMGU", 7, 24, true, -1);
    test_sites(sites, 8, "8mer_MMGU", 8, 24, true, -1);
    test_sites(sites, 9, "8mer_MMGU", 9, 24, true, -1);
    test_sites(sites, 10, "7mer_GUM", 10, 24, true, 1);
    test_sites(sites, 11, "7mer_GUM", 11, 24, true, 1);
    test_sites(sites, 12, "7mer_MMGU", 12, 24, true, -1);
    test_sites(sites, 13, "7mer_MMGU", 13, 24, true, -1);
    test_sites(sites, 14, "8mer_MMGU", 14, 24, true, -1);
    test_sites(sites, 15, "8mer_MMGU", 15, 24, true, -1);
    test_sites(sites, 16, "7mer_GUT", 16, 24, true, 4);
    test_sites(sites, 17, "7mer_GUT", 17, 24, true, 4);
    test_sites(sites, 18, "7mer_MMGU", 18, 24, true, -1);
    test_sites(sites, 19, "7mer_MMGU", 19, 24, true, -1);
    test_sites(sites, 20, "8mer_MMGU", 20, 24, true, -1);
    test_sites(sites, 21, "8mer_MMGU", 21, 24, true, -1);
    test_sites(sites, 22, "7mer_GUT", 22, 24, true, 5);
    test_sites(sites, 23, "7mer_GUT", 23, 24, true, 5);
    test_sites(sites, 24, "7mer_MMGU", 24, 24, true, -1);
    test_sites(sites, 25, "7mer_MMGU", 25, 24, true, -1);

    test_sites(sites, 26, "7mer_GU+", 0, 14, true, 0);
    test_sites(sites, 27, "7mer_GU+", 1, 14, true, 0);
    test_sites(sites, 28, "7mer_GU+", 2, 14, true, 0);
    test_sites(sites, 29, "7mer_GU+", 3, 14, true, 0);
    test_sites(sites, 30, "7mer_GU+", 4, 14, true, 0);
    test_sites(sites, 31, "7mer_GU+", 5, 14, true, 0);
    test_sites(sites, 32, "7mer_GU+", 6, 14, true, 0);
    test_sites(sites, 33, "7mer_GU+", 7, 14, true, 0);
    test_sites(sites, 34, "7mer_GU+", 8, 14, true, 0);
    test_sites(sites, 35, "7mer_GU+", 9, 14, true, 0);
    test_sites(sites, 36, "7mer_GU+", 10, 14, true, 0);
    test_sites(sites, 37, "7mer_GU+", 11, 14, true, 0);
    test_sites(sites, 38, "7mer_GU+", 12, 14, true, 0);
    test_sites(sites, 39, "7mer_GU+", 13, 14, true, 0);
    test_sites(sites, 40, "7mer_GU+", 14, 14, true, 0);
    test_sites(sites, 41, "7mer_GU+", 15, 14, true, 0);
    test_sites(sites, 42, "7mer_GU+", 16, 14, true, 0);
    test_sites(sites, 43, "7mer_GU+", 17, 14, true, 0);
    test_sites(sites, 44, "7mer_GU+", 18, 14, true, 0);
    test_sites(sites, 45, "7mer_GU+", 19, 14, true, 0);
    test_sites(sites, 46, "7mer_GU+", 20, 14, true, 0);
    test_sites(sites, 47, "7mer_GU+", 21, 14, true, 0);
    test_sites(sites, 48, "7mer_GU+", 22, 14, true, 0);
    test_sites(sites, 49, "7mer_GU+", 23, 14, true, 0);
    test_sites(sites, 50, "7mer_GU+", 24, 14, true, 0);
    test_sites(sites, 51, "7mer_GU+", 25, 14, true, 0);
    test_sites(sites, 52, "7mer_GU+", 26, 14, true, 0);
    test_sites(sites, 53, "7mer_GU+", 27, 14, true, 0);
    test_sites(sites, 54, "7mer_GU+", 28, 14, true, 0);
    test_sites(sites, 55, "7mer_GU+", 29, 14, true, 0);
    test_sites(sites, 56, "7mer_GU+", 30, 14, true, 0);
    test_sites(sites, 57, "7mer_GU+", 31, 14, true, 0);
    test_sites(sites, 58, "7mer_GU+", 32, 14, true, 0);
    test_sites(sites, 59, "7mer_GU+", 33, 14, true, 0);
    test_sites(sites, 60, "7mer_GU+", 34, 14, true, 0);
    test_sites(sites, 61, "7mer_GU+", 35, 14, true, 0);
    test_sites(sites, 62, "7mer_GU+", 36, 14, true, 0);
    test_sites(sites, 63, "7mer_GU+", 37, 14, true, 0);
    test_sites(sites, 64, "7mer_GU+", 38, 14, true, 0);
    test_sites(sites, 65, "7mer_GU+", 39, 14, true, 0);
    test_sites(sites, 66, "7mer_GU+", 40, 14, true, 0);
    test_sites(sites, 67, "7mer_GU+", 41, 14, true, 0);
    test_sites(sites, 68, "7mer_GU+", 42, 14, true, 0);
    test_sites(sites, 69, "7mer_GU+", 43, 14, true, 0);
    test_sites(sites, 70, "7mer_GU+", 44, 14, true, 0);
    test_sites(sites, 71, "7mer_GU+", 45, 14, true, 0);
    test_sites(sites, 72, "7mer_GU+", 46, 14, true, 0);
    test_sites(sites, 73, "7mer_GU+", 47, 14, true, 0);
    test_sites(sites, 74, "7mer_GU+", 48, 14, true, 0);
    test_sites(sites, 75, "7mer_GU+", 49, 14, true, 0);
    test_sites(sites, 76, "7mer_GU+", 50, 14, true, 0);
    test_sites(sites, 77, "7mer_GU+", 51, 14, true, 0);
    test_sites(sites, 78, "7mer_GU+", 52, 14, true, 0);
    test_sites(sites, 79, "7mer_GU+", 53, 14, true, 0);
    test_sites(sites, 80, "7mer_GU+", 54, 14, true, 0);
    test_sites(sites, 81, "7mer_GU+", 55, 14, true, 0);
    test_sites(sites, 82, "7mer_GU+", 56, 14, true, 0);
    test_sites(sites, 83, "7mer_GU+", 57, 14, true, 0);
    test_sites(sites, 84, "7mer_GU+", 58, 14, true, 0);
    test_sites(sites, 85, "7mer_GU+", 59, 14, true, 0);
    test_sites(sites, 86, "7mer_GU+", 60, 14, true, 0);
    test_sites(sites, 87, "7mer_GU+", 61, 14, true, 0);
    test_sites(sites, 88, "7mer_GU+", 62, 14, true, 0);
    test_sites(sites, 89, "7mer_GU+", 63, 14, true, 0);
    test_sites(sites, 90, "7mer_GU+", 64, 14, true, 0);
    test_sites(sites, 91, "7mer_GU+", 65, 14, true, 0);
    test_sites(sites, 92, "7mer_GU+", 66, 14, true, 0);
    test_sites(sites, 93, "7mer_GU+", 67, 14, true, 0);
    test_sites(sites, 94, "7mer_GU+", 68, 14, true, 0);
    test_sites(sites, 95, "7mer_GU+", 69, 14, true, 0);
    test_sites(sites, 96, "7mer_GU+", 70, 14, true, 0);
    test_sites(sites, 97, "7mer_GU+", 71, 14, true, 0);
    test_sites(sites, 98, "7mer_GU+", 72, 14, true, 0);
    test_sites(sites, 99, "7mer_GU+", 73, 14, true, 0);
    test_sites(sites, 100, "7mer_GU+", 74, 14, true, 0);
    test_sites(sites, 101, "7mer_GU+", 75, 14, true, 0);
    test_sites(sites, 102, "7mer_GU+", 76, 14, true, 0);
    test_sites(sites, 103, "7mer_GU+", 77, 14, true, 0);
    test_sites(sites, 104, "7mer_GU+", 78, 14, true, 0);
    test_sites(sites, 105, "7mer_GU+", 79, 14, true, 0);
    test_sites(sites, 106, "7mer_GU+", 80, 14, true, 0);
    test_sites(sites, 107, "7mer_GU+", 81, 14, true, 0);
    test_sites(sites, 108, "7mer_GU+", 82, 14, true, 0);
    test_sites(sites, 109, "7mer_GU+", 83, 14, true, 0);
    test_sites(sites, 110, "7mer_GU+", 84, 14, true, 0);
    test_sites(sites, 111, "7mer_GU+", 85, 14, true, 0);
    test_sites(sites, 112, "7mer_GU+", 86, 14, true, 0);
    test_sites(sites, 113, "7mer_GU+", 87, 14, true, 0);
    test_sites(sites, 114, "7mer_GU+", 88, 14, true, 0);
    test_sites(sites, 115, "7mer_GU+", 89, 14, true, 0);
    test_sites(sites, 116, "7mer_GU+", 90, 14, true, 0);
    test_sites(sites, 117, "7mer_GU+", 91, 14, true, 0);
    test_sites(sites, 118, "7mer_GU+", 92, 14, true, 0);
    test_sites(sites, 119, "7mer_GU+", 93, 14, true, 0);

    test_sites(sites, 120, "8mer_MMGU", 26, 24, true, 0);
    test_sites(sites, 121, "8mer_MMGU", 27, 24, true, 0);
    test_sites(sites, 122, "8mer_MMGU", 28, 24, true, 1);
    test_sites(sites, 123, "8mer_MMGU", 29, 24, true, 1);
    test_sites(sites, 124, "8mer_MMGU", 30, 24, true, 2);
    test_sites(sites, 125, "8mer_MMGU", 31, 24, true, 2);
    test_sites(sites, 126, "8mer_MMGU", 32, 24, true, 3);
    test_sites(sites, 127, "8mer_MMGU", 33, 24, true, 3);
    test_sites(sites, 128, "8mer_MMGU", 34, 24, true, 4);
    test_sites(sites, 129, "8mer_MMGU", 35, 24, true, 4);
    test_sites(sites, 130, "8mer_MMGU", 36, 24, true, 5);
    test_sites(sites, 131, "8mer_MMGU", 37, 24, true, 5);
    test_sites(sites, 132, "8mer_MMGU", 38, 24, true, 1);
    test_sites(sites, 133, "8mer_MMGU", 39, 24, true, 1);
    test_sites(sites, 134, "7mer_MMGU", 40, 24, true, 1);
    test_sites(sites, 135, "7mer_MMGU", 41, 24, true, 1);
    test_sites(sites, 136, "8mer_MMGU", 42, 24, true, 2);
    test_sites(sites, 137, "8mer_MMGU", 43, 24, true, 2);
    test_sites(sites, 138, "7mer_MMGU", 44, 24, true, 2);
    test_sites(sites, 139, "7mer_MMGU", 45, 24, true, 2);
    test_sites(sites, 140, "8mer_MMGU", 46, 24, true, 3);
    test_sites(sites, 141, "8mer_MMGU", 47, 24, true, 3);
    test_sites(sites, 142, "7mer_MMGU", 48, 24, true, 3);
    test_sites(sites, 143, "7mer_MMGU", 49, 24, true, 3);
    test_sites(sites, 144, "8mer_MMGU", 50, 24, true, 4);
    test_sites(sites, 145, "8mer_MMGU", 51, 24, true, 4);
    test_sites(sites, 146, "7mer_MMGU", 52, 24, true, 4);
    test_sites(sites, 147, "7mer_MMGU", 53, 24, true, 4);
    test_sites(sites, 148, "8mer_MMGU", 54, 24, true, 5);
    test_sites(sites, 149, "8mer_MMGU", 55, 24, true, 5);
    test_sites(sites, 150, "7mer_MMGU", 56, 24, true, 5);
    test_sites(sites, 151, "7mer_MMGU", 57, 24, true, 5);
    test_sites(sites, 152, "8mer_MMGU", 58, 24, true, 0);
    test_sites(sites, 153, "8mer_MMGU", 59, 24, true, 0);
    test_sites(sites, 154, "7mer_MMGU", 60, 24, true, 0);
    test_sites(sites, 155, "7mer_MMGU", 61, 24, true, 0);
    test_sites(sites, 156, "8mer_MMGU", 62, 24, true, 2);
    test_sites(sites, 157, "8mer_MMGU", 63, 24, true, 2);
    test_sites(sites, 158, "7mer_MMGU", 64, 24, true, 2);
    test_sites(sites, 159, "7mer_MMGU", 65, 24, true, 2);
    test_sites(sites, 160, "8mer_MMGU", 66, 24, true, 3);
    test_sites(sites, 161, "8mer_MMGU", 67, 24, true, 3);
    test_sites(sites, 162, "7mer_MMGU", 68, 24, true, 3);
    test_sites(sites, 163, "7mer_MMGU", 69, 24, true, 3);
    test_sites(sites, 164, "8mer_MMGU", 70, 24, true, 4);
    test_sites(sites, 165, "8mer_MMGU", 71, 24, true, 4);
    test_sites(sites, 166, "7mer_MMGU", 72, 24, true, 4);
    test_sites(sites, 167, "7mer_MMGU", 73, 24, true, 4);
    test_sites(sites, 168, "8mer_MMGU", 74, 24, true, 5);
    test_sites(sites, 169, "8mer_MMGU", 75, 24, true, 5);
    test_sites(sites, 170, "7mer_MMGU", 76, 24, true, 5);
    test_sites(sites, 171, "7mer_MMGU", 77, 24, true, 5);
    test_sites(sites, 172, "8mer_MMGU", 78, 24, true, 0);
    test_sites(sites, 173, "8mer_MMGU", 79, 24, true, 0);
    test_sites(sites, 174, "7mer_MMGU", 80, 24, true, 0);
    test_sites(sites, 175, "7mer_MMGU", 81, 24, true, 0);
    test_sites(sites, 176, "8mer_MMGU", 82, 24, true, 2);
    test_sites(sites, 177, "8mer_MMGU", 83, 24, true, 2);
    test_sites(sites, 178, "7mer_MMGU", 84, 24, true, 2);
    test_sites(sites, 179, "7mer_MMGU", 85, 24, true, 2);
    test_sites(sites, 180, "8mer_MMGU", 86, 24, true, 3);
    test_sites(sites, 181, "8mer_MMGU", 87, 24, true, 3);
    test_sites(sites, 182, "7mer_MMGU", 88, 24, true, 3);
    test_sites(sites, 183, "7mer_MMGU", 89, 24, true, 3);
    test_sites(sites, 184, "8mer_MMGU", 90, 24, true, 5);
    test_sites(sites, 185, "8mer_MMGU", 91, 24, true, 5);
    test_sites(sites, 186, "7mer_MMGU", 92, 24, true, 5);
    test_sites(sites, 187, "7mer_MMGU", 93, 24, true, 5);

    test_sites(sites, 188, "BT", 22, 23, false, 0);
    test_sites(sites, 189, "BT", 23, 23, false, 0);
    test_sites(sites, 190, "BT", 24, 23, false, 0);
    test_sites(sites, 191, "BT", 25, 23, false, 0);
    test_sites(sites, 192, "BT", 0, 23, false, 0);
    test_sites(sites, 193, "BT", 1, 23, false, 0);
    test_sites(sites, 194, "BT", 36, 23, false, 0);
    test_sites(sites, 195, "BT", 37, 23, false, 0);
    test_sites(sites, 196, "BT", 20, 23, false, 0);
    test_sites(sites, 197, "BT", 21, 23, false, 0);

    test_sites(sites, 198, "7mer_BT", 6, 23, true, 1);
    test_sites(sites, 199, "7mer_BT", 7, 23, true, 1);

    test_sites(sites, 200, "8mer_BT", 82, 24, true, 2);
    test_sites(sites, 201, "BT", 83, 24, false, 0);
    test_sites(sites, 202, "7mer_BT", 84, 24, true, 2);
    test_sites(sites, 203, "BT", 85, 24, false, 0);

    test_sites(sites, 204, "8mer_BT", 86, 24, true, 3);
    test_sites(sites, 205, "BT", 87, 24, false, 0);
    test_sites(sites, 206, "7mer_BT", 88, 24, true, 3);
    test_sites(sites, 207, "BT", 89, 24, false, 0);

    test_sites(sites, 208, "BT", 34, 24, false, 0);
    test_sites(sites, 209, "BT", 35, 24, false, 0);
    test_sites(sites, 210, "BT", 36, 24, false, 0);
    test_sites(sites, 211, "BT", 37, 24, false, 0);
}
}
