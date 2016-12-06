#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_main_io.hpp"
#include "test_site.hpp"

namespace {

    class Site05BT1 : public TestSiteMR3AS
    {
    protected:
        Site05BT1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_05_bt_1.fasta";
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

    TEST_F(Site05BT1, mir124_bt) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(2u, sites.get_length());

        test_sites(sites, 0, "7mer_BT", 0, 17, true, 0);
        test_sites(sites, 1, "BT", 0, 18, false, 0);
    }

}
