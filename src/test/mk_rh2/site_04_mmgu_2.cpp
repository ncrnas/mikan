#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

    class Site04MMGU2 : public TestSiteRH2
    {
    protected:
        Site04MMGU2() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_04_mmgu_2.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_rh2/";

            mSeedDef1 = "6mGU1";
            mOverlapDef = "orig";
        }

        typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TIndexQGram TIdx;
        typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TFinder TFin;
        typedef rh2mfe::RH2SeedSites<rh2mfe::TRNATYPE> TSit;

    };

    TEST_F(Site04MMGU2, mir1_mmgu) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef1, mOverlapDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(26u, sites.get_length());

        test_sites2(sites, 0, "6mer", 0, 24, true);
        test_sites2(sites, 1, "6mer", 1, 24, true);

        test_sites2(sites, 2, "6mer_GUT", 2, 24, true);
        test_sites2(sites, 3, "6mer_GUT", 3, 24, true);
        test_sites2(sites, 4, "6mer_GUT", 4, 24, true);
        test_sites2(sites, 5, "6mer_GUT", 5, 24, true);
        test_sites2(sites, 6, "6mer_GUT", 6, 24, true);
        test_sites2(sites, 7, "6mer_GUT", 7, 24, true);

        test_sites2(sites, 8, "6mer_GUM", 8, 24, true);
        test_sites2(sites, 9, "6mer_GUM", 9, 24, true);
        test_sites2(sites, 10, "6mer_GUM", 10, 24, true);
        test_sites2(sites, 11, "6mer_GUM", 11, 24, true);
        test_sites2(sites, 12, "6mer_GUM", 12, 24, true);
        test_sites2(sites, 13, "6mer_GUM", 13, 24, true);

        test_sites2(sites, 14, "6mer_GUT", 14, 24, true);
        test_sites2(sites, 15, "6mer_GUT", 15, 24, true);
        test_sites2(sites, 16, "6mer_GUT", 16, 24, true);
        test_sites2(sites, 17, "6mer_GUT", 17, 24, true);
        test_sites2(sites, 18, "6mer_GUT", 18, 24, true);
        test_sites2(sites, 19, "6mer_GUT", 19, 24, true);

        test_sites2(sites, 20, "6mer_GUT", 20, 24, true);
        test_sites2(sites, 21, "6mer_GUT", 21, 24, true);
        test_sites2(sites, 22, "6mer_GUT", 22, 24, true);
        test_sites2(sites, 23, "6mer_GUT", 23, 24, true);
        test_sites2(sites, 24, "6mer_GUT", 24, 24, true);
        test_sites2(sites, 25, "6mer_GUT", 25, 24, true);
    }

}
