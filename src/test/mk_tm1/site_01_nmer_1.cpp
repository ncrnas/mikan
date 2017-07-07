#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site01Nmer1 : public TestSiteTM1 {
protected:
    Site01Nmer1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tm1/";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef tm1p::TM1SeedSites<mikan::TRNATYPE> TSit;

};

TEST_F(Site01Nmer1, mir124) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[0]);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(55u, sites.get_length());

    test_sites2(sites, 0, "6mer", 0, 0, true);
    test_sites2(sites, 1, "6mer", 1, 1, true);
    test_sites2(sites, 2, "6mer", 2, 13, true);
    test_sites2(sites, 3, "6mer", 3, 14, true);
    test_sites2(sites, 4, "6mer", 4, 15, true);
    test_sites2(sites, 5, "6mer", 5, 18, true);
    test_sites2(sites, 6, "6mer", 6, 19, true);
    test_sites2(sites, 7, "6mer", 7, 20, true);
    test_sites2(sites, 8, "6mer", 8, 32, true);
    test_sites2(sites, 9, "6mer", 9, 33, true);
    test_sites2(sites, 10, "6mer", 10, 34, true);

    test_sites2(sites, 11, "7mer-A1", 11, 0, true);
    test_sites2(sites, 12, "7mer-A1", 12, 1, true);
    test_sites2(sites, 13, "7mer-A1", 13, 13, true);
    test_sites2(sites, 14, "7mer-A1", 14, 14, true);
    test_sites2(sites, 15, "7mer-A1", 15, 15, true);
    test_sites2(sites, 16, "7mer-A1", 16, 18, true);
    test_sites2(sites, 17, "7mer-A1", 17, 19, true);
    test_sites2(sites, 18, "7mer-A1", 18, 20, true);
    test_sites2(sites, 19, "7mer-A1", 19, 32, true);
    test_sites2(sites, 20, "7mer-A1", 20, 33, true);
    test_sites2(sites, 21, "6mer", 21, 34, true);

    test_sites2(sites, 22, "7mer-m8", 22, 1, true);
    test_sites2(sites, 23, "7mer-m8", 23, 13, true);
    test_sites2(sites, 24, "7mer-m8", 24, 14, true);
    test_sites2(sites, 25, "7mer-m8", 25, 15, true);
    test_sites2(sites, 26, "7mer-m8", 26, 18, true);
    test_sites2(sites, 27, "7mer-m8", 27, 19, true);
    test_sites2(sites, 28, "7mer-m8", 28, 20, true);
    test_sites2(sites, 29, "7mer-m8", 29, 21, true);
    test_sites2(sites, 30, "7mer-m8", 30, 32, true);
    test_sites2(sites, 31, "7mer-m8", 31, 33, true);
    test_sites2(sites, 32, "7mer-m8", 32, 34, true);

    test_sites2(sites, 33, "8mer", 33, 1, true);
    test_sites2(sites, 34, "8mer", 34, 13, true);
    test_sites2(sites, 35, "8mer", 35, 14, true);
    test_sites2(sites, 36, "8mer", 36, 15, true);
    test_sites2(sites, 37, "8mer", 37, 18, true);
    test_sites2(sites, 38, "8mer", 38, 19, true);
    test_sites2(sites, 39, "8mer", 39, 20, true);
    test_sites2(sites, 40, "8mer", 40, 21, true);
    test_sites2(sites, 41, "8mer", 41, 32, true);
    test_sites2(sites, 42, "8mer", 42, 33, true);
    test_sites2(sites, 43, "7mer-m8", 43, 34, true);

    test_sites2(sites, 44, "7mer-m8", 44, 2, true);
    test_sites2(sites, 45, "7mer-m8", 45, 13, true);
    test_sites2(sites, 46, "7mer-m8", 46, 14, true);
    test_sites2(sites, 47, "7mer-m8", 47, 15, true);
    test_sites2(sites, 48, "7mer-m8", 48, 18, true);
    test_sites2(sites, 49, "7mer-m8", 49, 19, true);
    test_sites2(sites, 50, "7mer-m8", 50, 20, true);
    test_sites2(sites, 51, "7mer-m8", 51, 21, true);
    test_sites2(sites, 52, "7mer-m8", 52, 32, true);
    test_sites2(sites, 53, "7mer-m8", 53, 33, true);
    test_sites2(sites, 54, "7mer-m8", 54, 34, true);
}

}
