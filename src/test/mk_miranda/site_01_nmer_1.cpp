#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

    class Site01Nmer1 : public TestSiteMR3AS
    {
    protected:
        Site01Nmer1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_01_nmer_1.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_miranda/";

            resize(mSeedDef, 6);
            mSeedDef[0] = 'Y';
            mSeedDef[1] = 'Y';
            mSeedDef[2] = 'Y';
            mSeedDef[3] = "0";
            mSeedDef[4] = "0:0";
            mSeedDef[5] = "0";
        }

        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TIndexQGram TIdx;
        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TFinder TFin;
        typedef mr3as::MR3SeedSites<mr3as::TRNATYPE> TSit;

    };

    TEST_F(Site01Nmer1, mir124_8mer) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(55u, sites.get_length());

        test_sites(sites, 0, "", 0, 0, false, 0);
        test_sites(sites, 1, "6mer", 1, 1, true, 0);
        test_sites(sites, 2, "6mer", 2, 13, true, 0);
        test_sites(sites, 3, "6mer", 3, 14, true, 0);
        test_sites(sites, 4, "6mer", 4, 15, true, 0);
        test_sites(sites, 5, "6mer", 5, 18, true, 0);
        test_sites(sites, 6, "6mer", 6, 19, true, 0);
        test_sites(sites, 7, "6mer", 7, 20, true, 0);
        test_sites(sites, 8, "6mer", 8, 32, true, 0);
        test_sites(sites, 9, "6mer", 9, 33, true, 0);
        test_sites(sites, 10, "", 10, 34, false, 0);

        test_sites(sites, 11, "", 11, 0, false, 0);
        test_sites(sites, 12, "6mer", 12, 1, true, 0);
        test_sites(sites, 13, "6mer", 13, 13, true, 0);
        test_sites(sites, 14, "6mer", 14, 14, true, 0);
        test_sites(sites, 15, "6mer", 15, 15, true, 0);
        test_sites(sites, 16, "6mer", 16, 18, true, 0);
        test_sites(sites, 17, "6mer", 17, 19, true, 0);
        test_sites(sites, 18, "6mer", 18, 20, true, 0);
        test_sites(sites, 19, "6mer", 19, 32, true, 0);
        test_sites(sites, 20, "6mer", 20, 33, true, 0);
        test_sites(sites, 21, "", 21, 34, false, 0);

        test_sites(sites, 22, "7mer", 22, 1, true, 0);
        test_sites(sites, 23, "7mer", 23, 13, true, 0);
        test_sites(sites, 24, "7mer", 24, 14, true, 0);
        test_sites(sites, 25, "7mer", 25, 15, true, 0);
        test_sites(sites, 26, "7mer", 26, 18, true, 0);
        test_sites(sites, 27, "7mer", 27, 19, true, 0);
        test_sites(sites, 28, "7mer", 28, 20, true, 0);
        test_sites(sites, 29, "7mer", 29, 21, true, 0);
        test_sites(sites, 30, "7mer", 30, 32, true, 0);
        test_sites(sites, 31, "7mer", 31, 33, true, 0);
        test_sites(sites, 32, "", 32, 34, false, 0);

        test_sites(sites, 33, "7mer", 33, 1, true, 0);
        test_sites(sites, 34, "7mer", 34, 13, true, 0);
        test_sites(sites, 35, "7mer", 35, 14, true, 0);
        test_sites(sites, 36, "7mer", 36, 15, true, 0);
        test_sites(sites, 37, "7mer", 37, 18, true, 0);
        test_sites(sites, 38, "7mer", 38, 19, true, 0);
        test_sites(sites, 39, "7mer", 39, 20, true, 0);
        test_sites(sites, 40, "7mer", 40, 21, true, 0);
        test_sites(sites, 41, "7mer", 41, 32, true, 0);
        test_sites(sites, 42, "7mer", 42, 33, true, 0);
        test_sites(sites, 43, "", 43, 34, false, 0);

        test_sites(sites, 44, "8mer", 44, 2, true, 0);
        test_sites(sites, 45, "8mer", 45, 13, true, 0);
        test_sites(sites, 46, "8mer", 46, 14, true, 0);
        test_sites(sites, 47, "8mer", 47, 15, true, 0);
        test_sites(sites, 48, "8mer", 48, 18, true, 0);
        test_sites(sites, 49, "8mer", 49, 19, true, 0);
        test_sites(sites, 50, "8mer", 50, 20, true, 0);
        test_sites(sites, 51, "8mer", 51, 21, true, 0);
        test_sites(sites, 52, "8mer", 52, 32, true, 0);
        test_sites(sites, 53, "8mer", 53, 33, true, 0);
        test_sites(sites, 54, "", 54, 34, false, 0);
    }

    TEST_F(Site01Nmer1, mir124_7mer) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        mSeedDef[2] = 'N';
        int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(55u, sites.get_length());

        test_sites(sites, 8, "6mer", 8, 32, true, 0);
        test_sites(sites, 19, "6mer", 19, 32, true, 0);
        test_sites(sites, 30, "7mer", 30, 32, true, 0);
        test_sites(sites, 41, "7mer", 41, 32, true, 0);
        test_sites(sites, 52, "7mer", 52, 32, true, 0);
    }

    TEST_F(Site01Nmer1, mir124_6mer) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        mSeedDef[1] = 'N';
        mSeedDef[2] = 'N';
        int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(55u, sites.get_length());

        test_sites(sites, 8, "6mer", 8, 32, true, 0);
        test_sites(sites, 19, "6mer", 19, 32, true, 0);
        test_sites(sites, 30, "6mer", 30, 32, true, 0);
        test_sites(sites, 41, "6mer", 41, 32, true, 0);
        test_sites(sites, 52, "6mer", 52, 32, true, 0);
    }

}
