#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

    class Site02GU1 : public TestSiteTM1
    {
    protected:
        Site02GU1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_02_gu_1.fasta";
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

    TEST_F(Site02GU1, mir124_gu) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[0]);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(11u, sites.get_length());

        test_sites2(sites, 0, "7mer-m8", 0, 24, true);

        test_sites2(sites, 1, "6mer", 1, 24, true);
        test_sites2(sites, 2, "7mer-m8", 2, 24, true);
        test_sites2(sites, 3, "8mer", 3, 24, true);
        test_sites2(sites, 4, "7mer-m8", 4, 24, true);
        test_sites2(sites, 5, "8mer", 5, 24, true);

        test_sites2(sites, 6, "6mer", 6, 24, true);
        test_sites2(sites, 7, "7mer-m8", 7, 24, true);
        test_sites2(sites, 8, "8mer", 8, 24, true);
        test_sites2(sites, 9, "7mer-m8", 9, 24, true);
        test_sites2(sites, 10, "8mer", 10, 24, true);
    }
}
