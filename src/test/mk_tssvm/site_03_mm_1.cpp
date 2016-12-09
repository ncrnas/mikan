#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

    class Site03MM1 : public TestSiteTSSVM
    {
    protected:
        Site03MM1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_03_mm_1.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_tssvm/";
        }

        typedef tssvm::TSSVMCore<tssvm::TRNATYPE>::TIndexQGram TIdx;
        typedef tssvm::TSSVMCore<tssvm::TRNATYPE>::TFinder TFin;
        typedef tssvm::TSSVMSeedSites<tssvm::TRNATYPE> TSit;

    };

    TEST_F(Site03MM1, mir124_mm) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[0]);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(16u, sites.get_length());

        test_sites(sites, 0, "6mer", 0, 24, true, 0);
        test_sites(sites, 1, "7mer-A1", 1, 24, true, 0);
        test_sites(sites, 2, "6mer", 2, 24, true, 0);
        test_sites(sites, 3, "7mer-A1", 3, 24, true, 0);

        test_sites(sites, 4, "LP", 37, 24, true, 1);
        test_sites(sites, 5, "LP", 39, 24, true, 1);
        test_sites(sites, 6, "LP", 31, 24, true, 2);
        test_sites(sites, 7, "LP", 33, 24, true, 2);
        test_sites(sites, 8, "LP", 25, 24, true, 3);
        test_sites(sites, 9, "LP", 27, 24, true, 3);
        test_sites(sites, 10, "LP", 19, 24, true, 4);
        test_sites(sites, 11, "LP", 21, 24, true, 4);
        test_sites(sites, 12, "LP", 13, 24, true, 5);
        test_sites(sites, 13, "LP", 15, 24, true, 5);
        test_sites(sites, 14, "LP", 7, 24, true, 6);
        test_sites(sites, 15, "LP", 9, 24, true, 6);
    }
}
