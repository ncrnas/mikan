#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

    class Site04MMGU1 : public TestSiteMR3AS
    {
    protected:
        Site04MMGU1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_04_mmgu_1.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_miranda/";

            resize(mSeedDef, 6);
            mSeedDef[0] = 'Y';
            mSeedDef[1] = 'Y';
            mSeedDef[2] = 'Y';
            mSeedDef[3] = "1";
            mSeedDef[4] = "1:1";
            mSeedDef[5] = "0";
        }

        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TIndexQGram TIdx;
        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TFinder TFin;
        typedef mr3as::MR3SeedSites<mr3as::TRNATYPE> TSit;

    };

    TEST_F(Site04MMGU1, mir124_mm8gu) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(66u, sites.get_length());

        test_sites(sites, 0, "7mer_MM", 0, 24, true, -1);
        test_sites(sites, 1, "7mer_MM", 1, 24, true, -1);

        test_sites(sites, 2, "8mer_MMGU", 2, 24, true, -1);
        test_sites(sites, 3, "8mer_MMGU", 3, 24, true, -1);
        test_sites(sites, 4, "7mer_GUT", 4, 24, true, 2);
        test_sites(sites, 5, "7mer_GUT", 5, 24, true, 2);
        test_sites(sites, 6, "7mer_MMGU", 6, 24, true, -1);
        test_sites(sites, 7, "7mer_MMGU", 7, 24, true, -1);

        test_sites(sites, 8, "8mer_MMGU", 8, 24, true, -1);
        test_sites(sites, 9, "8mer_MMGU", 9, 24, true, -1);
        test_sites(sites, 10, "7mer_GUT", 10, 24, true, 3);
        test_sites(sites, 11, "7mer_GUT", 11, 24, true, 3);
        test_sites(sites, 12, "7mer_MMGU", 12, 24, true, -1);
        test_sites(sites, 13, "7mer_MMGU", 13, 24, true, -1);

        test_sites(sites, 14, "7mer_MM", 14, 24, true, 0);
        test_sites(sites, 15, "7mer_MM", 15, 24, true, 0);
        test_sites(sites, 16, "7mer_MM", 16, 24, true, 1);
        test_sites(sites, 17, "7mer_MM", 17, 24, true, 1);
        test_sites(sites, 18, "7mer_MM", 18, 24, true, 2);
        test_sites(sites, 19, "7mer_MM", 19, 24, true, 2);
        test_sites(sites, 20, "7mer_MM", 20, 24, true, 3);
        test_sites(sites, 21, "7mer_MM", 21, 24, true, 3);
        test_sites(sites, 22, "7mer_MM", 22, 24, true, 4);
        test_sites(sites, 23, "7mer_MM", 23, 24, true, 4);
        test_sites(sites, 24, "7mer_MM", 24, 24, true, 5);
        test_sites(sites, 25, "7mer_MM", 25, 24, true, 5);

        test_sites(sites, 26, "8mer_MMGU", 26, 24, true, 0);
        test_sites(sites, 27, "8mer_MMGU", 27, 24, true, 0);
        test_sites(sites, 28, "7mer_MMGU", 28, 24, true, 0);
        test_sites(sites, 29, "7mer_MMGU", 29, 24, true, 0);
        test_sites(sites, 30, "8mer_MMGU", 30, 24, true, 1);
        test_sites(sites, 31, "8mer_MMGU", 31, 24, true, 1);
        test_sites(sites, 32, "7mer_MMGU", 32, 24, true, 1);
        test_sites(sites, 33, "7mer_MMGU", 33, 24, true, 1);
        test_sites(sites, 34, "8mer_MMGU", 34, 24, true, 3);
        test_sites(sites, 35, "8mer_MMGU", 35, 24, true, 3);
        test_sites(sites, 36, "7mer_MMGU", 36, 24, true, 3);
        test_sites(sites, 37, "7mer_MMGU", 37, 24, true, 3);
        test_sites(sites, 38, "8mer_MMGU", 38, 24, true, 4);
        test_sites(sites, 39, "8mer_MMGU", 39, 24, true, 4);
        test_sites(sites, 40, "7mer_MMGU", 40, 24, true, 4);
        test_sites(sites, 41, "7mer_MMGU", 41, 24, true, 4);
        test_sites(sites, 42, "8mer_MMGU", 42, 24, true, 5);
        test_sites(sites, 43, "8mer_MMGU", 43, 24, true, 5);
        test_sites(sites, 44, "7mer_MMGU", 44, 24, true, 5);
        test_sites(sites, 45, "7mer_MMGU", 45, 24, true, 5);

        test_sites(sites, 46, "8mer_MMGU", 46, 24, true, 0);
        test_sites(sites, 47, "8mer_MMGU", 47, 24, true, 0);
        test_sites(sites, 48, "7mer_MMGU", 48, 24, true, 0);
        test_sites(sites, 49, "7mer_MMGU", 49, 24, true, 0);
        test_sites(sites, 50, "8mer_MMGU", 50, 24, true, 1);
        test_sites(sites, 51, "8mer_MMGU", 51, 24, true, 1);
        test_sites(sites, 52, "7mer_MMGU", 52, 24, true, 1);
        test_sites(sites, 53, "7mer_MMGU", 53, 24, true, 1);
        test_sites(sites, 54, "8mer_MMGU", 54, 24, true, 2);
        test_sites(sites, 55, "8mer_MMGU", 55, 24, true, 2);
        test_sites(sites, 56, "7mer_MMGU", 56, 24, true, 2);
        test_sites(sites, 57, "7mer_MMGU", 57, 24, true, 2);
        test_sites(sites, 58, "8mer_MMGU", 58, 24, true, 4);
        test_sites(sites, 59, "8mer_MMGU", 59, 24, true, 4);
        test_sites(sites, 60, "7mer_MMGU", 60, 24, true, 4);
        test_sites(sites, 61, "7mer_MMGU", 61, 24, true, 4);
        test_sites(sites, 62, "8mer_MMGU", 62, 24, true, 5);
        test_sites(sites, 63, "8mer_MMGU", 63, 24, true, 5);
        test_sites(sites, 64, "7mer_MMGU", 64, 24, true, 5);
        test_sites(sites, 65, "7mer_MMGU", 65, 24, true, 5);
    }

    TEST_F(Site04MMGU1, mir124_mm7gu) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        mSeedDef[4] = "0:1";
        int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(66u, sites.get_length());

        test_sites(sites, 0, "8mer_MMGU", 0, 24, true, -1);
        test_sites(sites, 1, "8mer_MMGU", 1, 24, true, -1);

        test_sites(sites, 2, "8mer_MMGU", 2, 24, true, -1);
        test_sites(sites, 3, "8mer_MMGU", 3, 24, true, -1);
        test_sites(sites, 4, "7mer_GUT", 4, 24, true, 2);
        test_sites(sites, 5, "7mer_GUT", 5, 24, true, 2);
        test_sites(sites, 6, "GUT", 6, 24, false, 0);
        test_sites(sites, 7, "GUT", 7, 24, false, 0);

        test_sites(sites, 8, "8mer_MMGU", 8, 24, true, -1);
        test_sites(sites, 9, "8mer_MMGU", 9, 24, true, -1);
        test_sites(sites, 10, "7mer_GUT", 10, 24, true, 3);
        test_sites(sites, 11, "7mer_GUT", 11, 24, true, 3);
        test_sites(sites, 12, "GUT", 12, 24, false, 0);
        test_sites(sites, 13, "GUT", 13, 24, false, 0);

        test_sites(sites, 14, "8mer_MMGU", 14, 24, true, 0);
        test_sites(sites, 15, "8mer_MMGU", 15, 24, true, 0);
        test_sites(sites, 16, "8mer_MMGU", 16, 24, true, 1);
        test_sites(sites, 17, "8mer_MMGU", 17, 24, true, 1);
        test_sites(sites, 18, "8mer_MMGU", 18, 24, true, 2);
        test_sites(sites, 19, "8mer_MMGU", 19, 24, true, 2);
        test_sites(sites, 20, "8mer_MMGU", 20, 24, true, 3);
        test_sites(sites, 21, "8mer_MMGU", 21, 24, true, 3);
        test_sites(sites, 22, "8mer_MMGU", 22, 24, true, 4);
        test_sites(sites, 23, "8mer_MMGU", 23, 24, true, 4);
        test_sites(sites, 24, "8mer_MMGU", 24, 24, true, 5);
        test_sites(sites, 25, "8mer_MMGU", 25, 24, true, 5);

        test_sites(sites, 26, "8mer_MMGU", 26, 24, true, 0);
        test_sites(sites, 27, "8mer_MMGU", 27, 24, true, 0);
        test_sites(sites, 28, "MMGU", 28, 24, false, 0);
        test_sites(sites, 29, "MMGU", 29, 24, false, 0);
        test_sites(sites, 30, "8mer_MMGU", 30, 24, true, 1);
        test_sites(sites, 31, "8mer_MMGU", 31, 24, true, 1);
        test_sites(sites, 32, "MMGU", 32, 24, false, 0);
        test_sites(sites, 33, "MMGU", 33, 24, false, 0);
        test_sites(sites, 34, "8mer_MMGU", 34, 24, true, 3);
        test_sites(sites, 35, "8mer_MMGU", 35, 24, true, 3);
        test_sites(sites, 36, "MMGU", 36, 24, false, 0);
        test_sites(sites, 37, "MMGU", 37, 24, false, 0);
        test_sites(sites, 38, "8mer_MMGU", 38, 24, true, 4);
        test_sites(sites, 39, "8mer_MMGU", 39, 24, true, 4);
        test_sites(sites, 40, "MMGU", 40, 24, false, 0);
        test_sites(sites, 41, "MMGU", 41, 24, false, 0);
        test_sites(sites, 42, "8mer_MMGU", 42, 24, true, 5);
        test_sites(sites, 43, "8mer_MMGU", 43, 24, true, 5);
        test_sites(sites, 44, "MMGU", 44, 24, false, 0);
        test_sites(sites, 45, "MMGU", 45, 24, false, 0);

        test_sites(sites, 46, "8mer_MMGU", 46, 24, true, 0);
        test_sites(sites, 47, "8mer_MMGU", 47, 24, true, 0);
        test_sites(sites, 48, "MMGU", 48, 24, false, 0);
        test_sites(sites, 49, "MMGU", 49, 24, false, 0);
        test_sites(sites, 50, "8mer_MMGU", 50, 24, true, 1);
        test_sites(sites, 51, "8mer_MMGU", 51, 24, true, 1);
        test_sites(sites, 52, "MMGU", 52, 24, false, 0);
        test_sites(sites, 53, "MMGU", 53, 24, false, 0);
        test_sites(sites, 54, "8mer_MMGU", 54, 24, true, 2);
        test_sites(sites, 55, "8mer_MMGU", 55, 24, true, 2);
        test_sites(sites, 56, "MMGU", 56, 24, false, 0);
        test_sites(sites, 57, "MMGU", 57, 24, false, 0);
        test_sites(sites, 58, "8mer_MMGU", 58, 24, true, 4);
        test_sites(sites, 59, "8mer_MMGU", 59, 24, true, 4);
        test_sites(sites, 60, "MMGU", 60, 24, false, 0);
        test_sites(sites, 61, "MMGU", 61, 24, false, 0);
        test_sites(sites, 62, "8mer_MMGU", 62, 24, true, 5);
        test_sites(sites, 63, "8mer_MMGU", 63, 24, true, 5);
        test_sites(sites, 64, "MMGU", 64, 24, false, 0);
        test_sites(sites, 65, "MMGU", 65, 24, false, 0);
    }
}
