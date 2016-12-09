#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tssvm.hpp"

namespace {

    class Site01Nmer1 : public TestSiteTSSVM
    {
    protected:
        Site01Nmer1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_01_nmer_1.fasta";
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

    TEST_F(Site01Nmer1, mir124) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[0]);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(38u, sites.get_length());

        test_sites(sites, 0, "6mer", 4, 15, true, 0);
        test_sites(sites, 1, "6mer", 5, 18, true, 0);
        test_sites(sites, 2, "6mer", 6, 19, true, 0);
        test_sites(sites, 3, "6mer", 7, 20, true, 0);
        test_sites(sites, 4, "6mer", 8, 32, true, 0);
        test_sites(sites, 5, "6mer", 9, 33, true, 0);
        test_sites(sites, 6, "7mer-A1", 10, 34, true, 0);

        test_sites(sites, 7, "7mer-A1", 15, 15, true, 0);
        test_sites(sites, 8, "7mer-A1", 16, 18, true, 0);
        test_sites(sites, 9, "7mer-A1", 17, 19, true, 0);
        test_sites(sites, 10, "7mer-A1", 18, 20, true, 0);
        test_sites(sites, 11, "7mer-A1", 19, 32, true, 0);
        test_sites(sites, 12, "7mer-A1", 20, 33, true, 0);
        test_sites(sites, 13, "7mer-A1", 21, 34, true, 0);

        test_sites(sites, 14, "7mer-m8", 25, 15, true, 0);
        test_sites(sites, 15, "7mer-m8", 26, 18, true, 0);
        test_sites(sites, 16, "7mer-m8", 27, 19, true, 0);
        test_sites(sites, 17, "7mer-m8", 28, 20, true, 0);
        test_sites(sites, 18, "7mer-m8", 29, 21, true, 0);
        test_sites(sites, 19, "7mer-m8", 30, 32, true, 0);
        test_sites(sites, 20, "7mer-m8", 31, 33, true, 0);
        test_sites(sites, 21, "8mer", 32, 34, true, 0);

        test_sites(sites, 22, "8mer", 36, 15, true, 0);
        test_sites(sites, 23, "8mer", 37, 18, true, 0);
        test_sites(sites, 24, "8mer", 38, 19, true, 0);
        test_sites(sites, 25, "8mer", 39, 20, true, 0);
        test_sites(sites, 26, "8mer", 40, 21, true, 0);
        test_sites(sites, 27, "8mer", 41, 32, true, 0);
        test_sites(sites, 28, "8mer", 42, 33, true, 0);
        test_sites(sites, 29, "8mer", 43, 34, true, 0);

        test_sites(sites, 30, "7mer-m8", 47, 15, true, 0);
        test_sites(sites, 31, "7mer-m8", 48, 18, true, 0);
        test_sites(sites, 32, "7mer-m8", 49, 19, true, 0);
        test_sites(sites, 33, "7mer-m8", 50, 20, true, 0);
        test_sites(sites, 34, "7mer-m8", 51, 21, true, 0);
        test_sites(sites, 35, "7mer-m8", 52, 32, true, 0);
        test_sites(sites, 36, "7mer-m8", 53, 33, true, 0);
        test_sites(sites, 37, "8mer", 54, 34, true, 0);
    }

}
