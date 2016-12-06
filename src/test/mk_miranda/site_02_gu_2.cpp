#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

    class Site02GU2 : public TestSiteMR3AS
    {
    protected:
        Site02GU2() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_02_gu_2.fasta";
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
            mSeedDef[4] = "0:0";
            mSeedDef[5] = "0";
        }

        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TIndexQGram TIdx;
        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TFinder TFin;
        typedef mr3as::MR3SeedSites<mr3as::TRNATYPE> TSit;

    };

    TEST_F(Site02GU2, mir1_gu) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(40u, sites.get_length());

        test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
        test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);
        test_sites(sites, 2, "8mer_GUM", 2, 24, true, -1);
        test_sites(sites, 3, "8mer_GUM", 3, 24, true, -1);

        test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
        test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
        test_sites(sites, 6, "8mer_GUT", 6, 24, true, 0);
        test_sites(sites, 7, "8mer_GUT", 7, 24, true, 0);

        test_sites(sites, 8, "GUT", 8, 24, false, 0);
        test_sites(sites, 9, "GUT", 9, 24, false, 0);
        test_sites(sites, 10, "GUT", 10, 24, false, 0);
        test_sites(sites, 11, "GUT", 11, 24, false, 0);
        test_sites(sites, 12, "GUT", 12, 24, false, 0);

        test_sites(sites, 13, "7mer_GUM", 13, 24, true, 1);
        test_sites(sites, 14, "7mer_GUM", 14, 24, true, 1);
        test_sites(sites, 15, "8mer_GUM", 15, 24, true, 1);
        test_sites(sites, 16, "8mer_GUM", 16, 24, true, 1);

        test_sites(sites, 17, "GUM", 17, 24, false, 0);
        test_sites(sites, 18, "GUM", 18, 24, false, 0);
        test_sites(sites, 19, "GUM", 19, 24, false, 0);
        test_sites(sites, 20, "GUM", 20, 24, false, 0);
        test_sites(sites, 21, "GUM", 21, 24, false, 0);

        test_sites(sites, 22, "7mer_GUT", 22, 24, true, 4);
        test_sites(sites, 23, "7mer_GUT", 23, 24, true, 4);
        test_sites(sites, 24, "8mer_GUT", 24, 24, true, 4);
        test_sites(sites, 25, "8mer_GUT", 25, 24, true, 4);

        test_sites(sites, 26, "GUT", 26, 24, false, 0);
        test_sites(sites, 27, "GUT", 27, 24, false, 0);
        test_sites(sites, 28, "GUT", 28, 24, false, 0);
        test_sites(sites, 29, "GUT", 29, 24, false, 0);
        test_sites(sites, 30, "GUT", 30, 24, false, 0);

        test_sites(sites, 31, "7mer_GUT", 31, 24, true, 5);
        test_sites(sites, 32, "7mer_GUT", 32, 24, true, 5);
        test_sites(sites, 33, "8mer_GUT", 33, 24, true, 5);
        test_sites(sites, 34, "8mer_GUT", 34, 24, true, 5);

        test_sites(sites, 35, "GUT", 35, 24, false, 0);
        test_sites(sites, 36, "GUT", 36, 24, false, 0);
        test_sites(sites, 37, "GUT", 37, 24, false, 0);
        test_sites(sites, 38, "GUT", 38, 24, false, 0);
        test_sites(sites, 39, "GUT", 39, 24, false, 0);
    }

    TEST_F(Site02GU2, mir1_gu_plus) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        mSeedDef[3] = "+";
        int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(40u, sites.get_length());

        test_sites(sites, 0, "7mer_GUM", 0, 24, true, -1);
        test_sites(sites, 1, "7mer_GUM", 1, 24, true, -1);
        test_sites(sites, 2, "8mer_GUM", 2, 24, true, -1);
        test_sites(sites, 3, "8mer_GUM", 3, 24, true, -1);

        test_sites(sites, 4, "7mer_GUT", 4, 24, true, 0);
        test_sites(sites, 5, "7mer_GUT", 5, 24, true, 0);
        test_sites(sites, 6, "8mer_GUT", 6, 24, true, 0);
        test_sites(sites, 7, "8mer_GUT", 7, 24, true, 0);

        test_sites(sites, 8, "GUT", 8, 24, false, 0);
        test_sites(sites, 9, "7mer_GU+", 9, 24, true, 0);
        test_sites(sites, 10, "7mer_GU+", 10, 24, true, 0);
        test_sites(sites, 11, "8mer_GU+", 11, 24, true, 0);
        test_sites(sites, 12, "8mer_GU+", 12, 24, true, 0);

        test_sites(sites, 13, "7mer_GUM", 13, 24, true, 1);
        test_sites(sites, 14, "7mer_GUM", 14, 24, true, 1);
        test_sites(sites, 15, "8mer_GUM", 15, 24, true, 1);
        test_sites(sites, 16, "8mer_GUM", 16, 24, true, 1);

        test_sites(sites, 17, "GUM", 17, 24, false, 0);
        test_sites(sites, 18, "7mer_GU+", 18, 24, true, 0);
        test_sites(sites, 19, "7mer_GU+", 19, 24, true, 0);
        test_sites(sites, 20, "8mer_GU+", 20, 24, true, 0);
        test_sites(sites, 21, "8mer_GU+", 21, 24, true, 0);

        test_sites(sites, 22, "7mer_GUT", 22, 24, true, 4);
        test_sites(sites, 23, "7mer_GUT", 23, 24, true, 4);
        test_sites(sites, 24, "8mer_GUT", 24, 24, true, 4);
        test_sites(sites, 25, "8mer_GUT", 25, 24, true, 4);

        test_sites(sites, 26, "GUT", 26, 24, false, 0);
        test_sites(sites, 27, "7mer_GU+", 27, 24, true, 0);
        test_sites(sites, 28, "7mer_GU+", 28, 24, true, 0);
        test_sites(sites, 29, "8mer_GU+", 29, 24, true, 0);
        test_sites(sites, 30, "8mer_GU+", 30, 24, true, 0);

        test_sites(sites, 31, "7mer_GUT", 31, 24, true, 5);
        test_sites(sites, 32, "7mer_GUT", 32, 24, true, 5);
        test_sites(sites, 33, "8mer_GUT", 33, 24, true, 5);
        test_sites(sites, 34, "8mer_GUT", 34, 24, true, 5);

        test_sites(sites, 35, "GUT", 35, 24, false, 0);
        test_sites(sites, 36, "7mer_GU+", 36, 24, true, 0);
        test_sites(sites, 37, "7mer_GU+", 37, 24, true, 0);
        test_sites(sites, 38, "8mer_GU+", 38, 24, true, 0);
        test_sites(sites, 39, "8mer_GU+", 39, 24, true, 0);
    }
}
