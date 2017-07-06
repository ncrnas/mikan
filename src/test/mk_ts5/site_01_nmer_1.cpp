#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class Site01Nmer1 : public TestSiteTS5 {
protected:
    Site01Nmer1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_ts5/";
    }

    typedef ts5cs::TS5Core<ts5cs::TRNATYPE>::TIndexQGram TIdx;
    typedef ts5cs::TS5Core<ts5cs::TRNATYPE>::TFinder TFin;
    typedef ts5cs::TS5SeedSites<ts5cs::TRNATYPE> TSit;

};

TEST_F(Site01Nmer1, mir124) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder);

    int ret_val = sites.find_seed_sites(mirna_seqs[0]);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(55u, sites.get_length());

    test_sites3(sites, 0, 0, 0);
    test_sites3(sites, 1, 1, 1);
    test_sites3(sites, 2, 2, 13);
    test_sites3(sites, 3, 3, 14);
    test_sites3(sites, 4, 4, 15);
    test_sites3(sites, 5, 5, 18);
    test_sites3(sites, 6, 6, 19);
    test_sites3(sites, 7, 7, 20);
    test_sites3(sites, 8, 8, 32);
    test_sites3(sites, 9, 9, 33);
    test_sites3(sites, 10, 10, 34);

    test_sites3(sites, 11, 11, 0);
    test_sites3(sites, 12, 12, 1);
    test_sites3(sites, 13, 13, 13);
    test_sites3(sites, 14, 14, 14);
    test_sites3(sites, 15, 15, 15);
    test_sites3(sites, 16, 16, 18);
    test_sites3(sites, 17, 17, 19);
    test_sites3(sites, 18, 18, 20);
    test_sites3(sites, 19, 19, 32);
    test_sites3(sites, 20, 20, 33);
    test_sites3(sites, 21, 21, 34);

    test_sites3(sites, 22, 22, 1);
    test_sites3(sites, 23, 23, 13);
    test_sites3(sites, 24, 24, 14);
    test_sites3(sites, 25, 25, 15);
    test_sites3(sites, 26, 26, 18);
    test_sites3(sites, 27, 27, 19);
    test_sites3(sites, 28, 28, 20);
    test_sites3(sites, 29, 29, 21);
    test_sites3(sites, 30, 30, 32);
    test_sites3(sites, 31, 31, 33);
    test_sites3(sites, 32, 32, 34);

    test_sites3(sites, 33, 33, 1);
    test_sites3(sites, 34, 34, 13);
    test_sites3(sites, 35, 35, 14);
    test_sites3(sites, 36, 36, 15);
    test_sites3(sites, 37, 37, 18);
    test_sites3(sites, 38, 38, 19);
    test_sites3(sites, 39, 39, 20);
    test_sites3(sites, 40, 40, 21);
    test_sites3(sites, 41, 41, 32);
    test_sites3(sites, 42, 42, 33);
    test_sites3(sites, 43, 43, 34);

    test_sites3(sites, 44, 44, 2);
    test_sites3(sites, 45, 45, 13);
    test_sites3(sites, 46, 46, 14);
    test_sites3(sites, 47, 47, 15);
    test_sites3(sites, 48, 48, 18);
    test_sites3(sites, 49, 49, 19);
    test_sites3(sites, 50, 50, 20);
    test_sites3(sites, 51, 51, 21);
    test_sites3(sites, 52, 52, 32);
    test_sites3(sites, 53, 53, 33);
    test_sites3(sites, 54, 54, 34);
}

}
