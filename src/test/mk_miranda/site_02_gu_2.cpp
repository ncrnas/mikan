#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_main_io.hpp"
#include "test_site.hpp"

namespace {

    class Site02GU1 : public TestSiteMR3AS
    {
    protected:
        Site02GU1() {
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

    TEST_F(Site02GU1, mir2) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
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

        test_sites(sites, 8, "7mer_GUM", 8, 24, true, 1);
        test_sites(sites, 9, "7mer_GUM", 9, 24, true, 1);
        test_sites(sites, 10, "8mer_GUM", 10, 24, true, 1);
        test_sites(sites, 11, "8mer_GUM", 11, 24, true, 1);

        test_sites(sites, 12, "7mer_GUT", 12, 24, true, 4);
        test_sites(sites, 13, "7mer_GUT", 13, 24, true, 4);
        test_sites(sites, 14, "8mer_GUT", 14, 24, true, 4);
        test_sites(sites, 15, "8mer_GUT", 15, 24, true, 4);

        test_sites(sites, 16, "7mer_GUT", 16, 24, true, 5);
        test_sites(sites, 17, "7mer_GUT", 17, 24, true, 5);
        test_sites(sites, 18, "8mer_GUT", 18, 24, true, 5);
        test_sites(sites, 19, "8mer_GUT", 19, 24, true, 5);
    }
}
