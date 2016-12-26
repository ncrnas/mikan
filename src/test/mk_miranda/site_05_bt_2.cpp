#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

    class Site05BT2 : public TestSiteMR3AS
    {
    protected:
        Site05BT2() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_05_bt_2.fasta";
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
            mSeedDef[5] = "1";
        }

        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TIndexQGram TIdx;
        typedef mr3as::MR3Core<mr3as::TRNATYPE>::TFinder TFin;
        typedef mr3as::MR3SeedSites<mr3as::TRNATYPE> TSit;

    };

    TEST_F(Site05BT2, mir1_bt) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(25u, sites.get_length());

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

}
