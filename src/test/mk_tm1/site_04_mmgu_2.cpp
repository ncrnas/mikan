#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

    class Site04MMGU2 : public TestSiteTM1
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
            OMPATH = (char *)"mk_tm1/";
        }

        typedef tm1p::TM1Core<tm1p::TRNATYPE>::TIndexQGram TIdx;
        typedef tm1p::TM1Core<tm1p::TRNATYPE>::TFinder TFin;
        typedef tm1p::TM1SeedSites<tm1p::TRNATYPE> TSit;

    };

    TEST_F(Site04MMGU2, mir1_mmgu) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[1]);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(40u, sites.get_length());

        test_sites2(sites, 0, "6mer", 0, 24, true);
        test_sites2(sites, 1, "7mer-A1", 1, 24, true);

        test_sites2(sites, 2, "6mer", 36, 24, true);
        test_sites2(sites, 3, "6mer", 37, 24, true);

        test_sites2(sites, 4, "6mer", 20, 24, true);
        test_sites2(sites, 5, "7mer-A1", 21, 24, true);
        test_sites2(sites, 6, "7mer-m8", 22, 24, true);
        test_sites2(sites, 7, "8mer", 23, 24, true);
        test_sites2(sites, 8, "6mer", 24, 24, true);
        test_sites2(sites, 9, "7mer-A1", 25, 24, true);

        test_sites2(sites, 10, "6mer", 2, 24, true);
        test_sites2(sites, 11, "7mer-A1", 3, 24, true);
        test_sites2(sites, 12, "7mer-m8", 4, 24, true);
        test_sites2(sites, 13, "8mer", 5, 24, true);
        test_sites2(sites, 14, "6mer", 6, 24, true);
        test_sites2(sites, 15, "7mer-A1", 7, 24, true);

        test_sites2(sites, 16, "6mer", 54, 24, true);
        test_sites2(sites, 17, "6mer", 55, 24, true);
        test_sites2(sites, 18, "6mer", 56, 24, true);
        test_sites2(sites, 19, "6mer", 57, 24, true);

        test_sites2(sites, 20, "6mer", 8, 24, true);
        test_sites2(sites, 21, "7mer-A1", 9, 24, true);
        test_sites2(sites, 22, "7mer-m8", 10, 24, true);
        test_sites2(sites, 23, "8mer", 11, 24, true);
        test_sites2(sites, 24, "6mer", 12, 24, true);
        test_sites2(sites, 25, "7mer-A1", 13, 24, true);

        test_sites2(sites, 26, "6mer", 74, 24, true);
        test_sites2(sites, 27, "6mer", 75, 24, true);
        test_sites2(sites, 28, "6mer", 76, 24, true);
        test_sites2(sites, 29, "6mer", 77, 24, true);

        test_sites2(sites, 30, "6mer", 14, 24, true);
        test_sites2(sites, 31, "7mer-A1", 15, 24, true);
        test_sites2(sites, 32, "7mer-m8", 16, 24, true);
        test_sites2(sites, 33, "8mer", 17, 24, true);
        test_sites2(sites, 34, "6mer", 18, 24, true);
        test_sites2(sites, 35, "7mer-A1",19, 24, true);

        test_sites2(sites, 36, "6mer", 90, 24, true);
        test_sites2(sites, 37, "6mer", 91, 24, true);
        test_sites2(sites, 38, "6mer", 92, 24, true);
        test_sites2(sites, 39, "6mer", 93, 24, true);
    }

}
