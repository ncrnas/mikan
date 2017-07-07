#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site02GU2 : public TestSiteTM1 {
protected:
    Site02GU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tm1/";
    }

    typedef tm1p::TM1Core<mikan::TRNATYPE>::TIndexQGram TIdx;
    typedef tm1p::TM1Core<mikan::TRNATYPE>::TFinder TFin;
    typedef tm1p::TM1SeedSites<mikan::TRNATYPE> TSit;

};

TEST_F(Site02GU2, mir1_gu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[1]);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(40u, sites.get_length());

    test_sites2(sites, 0, "6mer", 0, 24, true);
    test_sites2(sites, 1, "7mer-A1", 1, 24, true);
    test_sites2(sites, 2, "6mer", 2, 24, true);
    test_sites2(sites, 3, "7mer-A1", 3, 24, true);

    test_sites2(sites, 4, "7mer-m8", 31, 24, true);
    test_sites2(sites, 5, "8mer", 32, 24, true);
    test_sites2(sites, 6, "7mer-m8", 33, 24, true);
    test_sites2(sites, 7, "8mer", 34, 24, true);

    test_sites2(sites, 8, "6mer", 35, 24, true);
    test_sites2(sites, 9, "6mer", 36, 24, true);
    test_sites2(sites, 10, "7mer-A1", 37, 24, true);
    test_sites2(sites, 11, "6mer", 38, 24, true);
    test_sites2(sites, 12, "7mer-A1", 39, 24, true);

    test_sites2(sites, 13, "7mer-m8", 4, 24, true);
    test_sites2(sites, 14, "8mer", 5, 24, true);
    test_sites2(sites, 15, "7mer-m8", 6, 24, true);
    test_sites2(sites, 16, "8mer", 7, 24, true);

    test_sites2(sites, 17, "6mer", 8, 24, true);
    test_sites2(sites, 18, "6mer", 9, 24, true);
    test_sites2(sites, 19, "7mer-A1", 10, 24, true);
    test_sites2(sites, 20, "6mer", 11, 24, true);
    test_sites2(sites, 21, "7mer-A1", 12, 24, true);

    test_sites2(sites, 22, "7mer-m8", 13, 24, true);
    test_sites2(sites, 23, "8mer", 14, 24, true);
    test_sites2(sites, 24, "7mer-m8", 15, 24, true);
    test_sites2(sites, 25, "8mer", 16, 24, true);

    test_sites2(sites, 26, "6mer", 17, 24, true);
    test_sites2(sites, 27, "6mer", 18, 24, true);
    test_sites2(sites, 28, "7mer-A1", 19, 24, true);
    test_sites2(sites, 29, "6mer", 20, 24, true);
    test_sites2(sites, 30, "7mer-A1", 21, 24, true);

    test_sites2(sites, 31, "7mer-m8", 22, 24, true);
    test_sites2(sites, 32, "8mer", 23, 24, true);
    test_sites2(sites, 33, "7mer-m8", 24, 24, true);
    test_sites2(sites, 34, "8mer", 25, 24, true);

    test_sites2(sites, 35, "6mer", 26, 24, true);
    test_sites2(sites, 36, "6mer", 27, 24, true);
    test_sites2(sites, 37, "7mer-A1", 28, 24, true);
    test_sites2(sites, 38, "6mer", 29, 24, true);
    test_sites2(sites, 39, "7mer-A1", 30, 24, true);
}
}
