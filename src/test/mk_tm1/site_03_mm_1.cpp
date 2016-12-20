#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

    class Site03MM1 : public TestSiteTM1
    {
    protected:
        Site03MM1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_03_mm_1.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_tm1/";
        }

        typedef tm1p::TM1Core<tm1p::TRNATYPE>::TIndexQGram TIdx;
        typedef tm1p::TM1Core<tm1p::TRNATYPE>::TFinder TFin;
        typedef tm1p::TM1SeedSites<tm1p::TRNATYPE> TSit;

    };

    TEST_F(Site03MM1, mir124_mm) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[0]);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(10u, sites.get_length());

        test_sites2(sites, 0, "6mer", 0, 24, true);
        test_sites2(sites, 1, "7mer-A1", 1, 24, true);
        test_sites2(sites, 2, "6mer", 2, 24, true);
        test_sites2(sites, 3, "7mer-A1", 3, 24, true);

        test_sites2(sites, 4, "", 34, 24, false);
        test_sites2(sites, 5, "", 35, 24, false);
        test_sites2(sites, 6, "6mer", 36, 24, true);
        test_sites2(sites, 7, "6mer", 37, 24, true);
        test_sites2(sites, 8, "6mer", 38, 24, true);
        test_sites2(sites, 9, "6mer", 39, 24, true);

    }
}