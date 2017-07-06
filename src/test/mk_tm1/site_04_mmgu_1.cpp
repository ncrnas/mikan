#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class Site04MMGU1 : public TestSiteTM1 {
protected:
    Site04MMGU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tm1/";
    }

    typedef tm1p::TM1Core<tm1p::TRNATYPE>::TIndexQGram TIdx;
    typedef tm1p::TM1Core<tm1p::TRNATYPE>::TFinder TFin;
    typedef tm1p::TM1SeedSites<tm1p::TRNATYPE> TSit;

};

TEST_F(Site04MMGU1, mir124_mmgu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[0]);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(24u, sites.get_length());

    test_sites2(sites, 0, "6mer", 0, 24, true);
    test_sites2(sites, 1, "7mer-A1", 1, 24, true);

    test_sites2(sites, 2, "6mer", 24, 24, true);
    test_sites2(sites, 3, "6mer", 25, 24, true);

    test_sites2(sites, 4, "6mer", 2, 24, true);
    test_sites2(sites, 5, "7mer-A1", 3, 24, true);
    test_sites2(sites, 6, "7mer-m8", 4, 24, true);
    test_sites2(sites, 7, "8mer", 5, 24, true);
    test_sites2(sites, 8, "6mer", 6, 24, true);
    test_sites2(sites, 9, "7mer-A1", 7, 24, true);

    test_sites2(sites, 10, "6mer", 42, 24, true);
    test_sites2(sites, 11, "6mer", 43, 24, true);
    test_sites2(sites, 12, "6mer", 44, 24, true);
    test_sites2(sites, 13, "6mer", 45, 24, true);

    test_sites2(sites, 14, "6mer", 8, 24, true);
    test_sites2(sites, 15, "7mer-A1", 9, 24, true);
    test_sites2(sites, 16, "7mer-m8", 10, 24, true);
    test_sites2(sites, 17, "8mer", 11, 24, true);
    test_sites2(sites, 18, "6mer", 12, 24, true);
    test_sites2(sites, 19, "7mer-A1", 13, 24, true);

    test_sites2(sites, 20, "6mer", 62, 24, true);
    test_sites2(sites, 21, "6mer", 63, 24, true);
    test_sites2(sites, 22, "6mer", 64, 24, true);
    test_sites2(sites, 23, "6mer", 65, 24, true);
}
}
