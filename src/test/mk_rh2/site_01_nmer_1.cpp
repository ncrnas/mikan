#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site01Nmer1 : public TestSiteRH2 {
protected:
    Site01Nmer1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        mSeedDef1 = "6mer";
        mOverlapDef = "orig";
    }

    typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TIndexQGram TIdx;
    typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TFinder TFin;
    typedef rh2mfe::RH2SeedSites<rh2mfe::TRNATYPE> TSit;

};

TEST_F(Site01Nmer1, mir124_6mer) {
    read_files(false);
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef1, mOverlapDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(55u, sites.get_length());

    test_sites2(sites, 0, "", 0, 0, false);
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

    test_sites2(sites, 11, "", 11, 0, false);
    test_sites2(sites, 12, "6mer", 12, 1, true);
    test_sites2(sites, 13, "6mer", 13, 13, true);
    test_sites2(sites, 14, "6mer", 14, 14, true);
    test_sites2(sites, 15, "6mer", 15, 15, true);
    test_sites2(sites, 16, "6mer", 16, 18, true);
    test_sites2(sites, 17, "6mer", 17, 19, true);
    test_sites2(sites, 18, "6mer", 18, 20, true);
    test_sites2(sites, 19, "6mer", 19, 32, true);
    test_sites2(sites, 20, "6mer", 20, 33, true);
    test_sites2(sites, 21, "6mer", 21, 34, true);

    test_sites2(sites, 22, "6mer", 22, 1, true);
    test_sites2(sites, 23, "6mer", 23, 13, true);
    test_sites2(sites, 24, "6mer", 24, 14, true);
    test_sites2(sites, 25, "6mer", 25, 15, true);
    test_sites2(sites, 26, "6mer", 26, 18, true);
    test_sites2(sites, 27, "6mer", 27, 19, true);
    test_sites2(sites, 28, "6mer", 28, 20, true);
    test_sites2(sites, 29, "6mer", 29, 21, true);
    test_sites2(sites, 30, "6mer", 30, 32, true);
    test_sites2(sites, 31, "6mer", 31, 33, true);
    test_sites2(sites, 32, "6mer", 32, 34, true);

    test_sites2(sites, 33, "6mer", 33, 1, true);
    test_sites2(sites, 34, "6mer", 34, 13, true);
    test_sites2(sites, 35, "6mer", 35, 14, true);
    test_sites2(sites, 36, "6mer", 36, 15, true);
    test_sites2(sites, 37, "6mer", 37, 18, true);
    test_sites2(sites, 38, "6mer", 38, 19, true);
    test_sites2(sites, 39, "6mer", 39, 20, true);
    test_sites2(sites, 40, "6mer", 40, 21, true);
    test_sites2(sites, 41, "6mer", 41, 32, true);
    test_sites2(sites, 42, "6mer", 42, 33, true);
    test_sites2(sites, 43, "6mer", 43, 34, true);

    test_sites2(sites, 44, "6mer", 44, 2, true);
    test_sites2(sites, 45, "6mer", 45, 13, true);
    test_sites2(sites, 46, "6mer", 46, 14, true);
    test_sites2(sites, 47, "6mer", 47, 15, true);
    test_sites2(sites, 48, "6mer", 48, 18, true);
    test_sites2(sites, 49, "6mer", 49, 19, true);
    test_sites2(sites, 50, "6mer", 50, 20, true);
    test_sites2(sites, 51, "6mer", 51, 21, true);
    test_sites2(sites, 52, "6mer", 52, 32, true);
    test_sites2(sites, 53, "6mer", 53, 33, true);
    test_sites2(sites, 54, "6mer", 54, 34, true);
}

TEST_F(Site01Nmer1, mir124_7mer) {
    read_files(false);
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef1 = "7mer";
    int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef1, mOverlapDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(55u, sites.get_length());

    test_sites2(sites, 0, "", 0, 0, false);
    test_sites2(sites, 1, "7mer", 1, 1, false);
    test_sites2(sites, 2, "7mer", 2, 13, false);
    test_sites2(sites, 3, "7mer", 3, 14, false);
    test_sites2(sites, 4, "7mer", 4, 15, false);
    test_sites2(sites, 5, "7mer", 5, 18, false);
    test_sites2(sites, 6, "7mer", 6, 19, false);
    test_sites2(sites, 7, "7mer", 7, 20, false);
    test_sites2(sites, 8, "7mer", 8, 32, false);
    test_sites2(sites, 9, "7mer", 9, 33, false);
    test_sites2(sites, 10, "7mer", 10, 34, false);

    test_sites2(sites, 11, "", 11, 0, false);
    test_sites2(sites, 12, "7mer", 12, 1, false);
    test_sites2(sites, 13, "7mer", 13, 13, false);
    test_sites2(sites, 14, "7mer", 14, 14, false);
    test_sites2(sites, 15, "7mer", 15, 15, false);
    test_sites2(sites, 16, "7mer", 16, 18, false);
    test_sites2(sites, 17, "7mer", 17, 19, false);
    test_sites2(sites, 18, "7mer", 18, 20, false);
    test_sites2(sites, 19, "7mer", 19, 32, false);
    test_sites2(sites, 20, "7mer", 20, 33, false);
    test_sites2(sites, 21, "7mer", 21, 34, false);

    test_sites2(sites, 22, "7mer", 22, 1, true);
    test_sites2(sites, 23, "7mer", 23, 13, true);
    test_sites2(sites, 24, "7mer", 24, 14, true);
    test_sites2(sites, 25, "7mer", 25, 15, true);
    test_sites2(sites, 26, "7mer", 26, 18, true);
    test_sites2(sites, 27, "7mer", 27, 19, true);
    test_sites2(sites, 28, "7mer", 28, 20, true);
    test_sites2(sites, 29, "7mer", 29, 21, true);
    test_sites2(sites, 30, "7mer", 30, 32, true);
    test_sites2(sites, 31, "7mer", 31, 33, true);
    test_sites2(sites, 32, "7mer", 32, 34, true);

    test_sites2(sites, 33, "7mer", 33, 1, true);
    test_sites2(sites, 34, "7mer", 34, 13, true);
    test_sites2(sites, 35, "7mer", 35, 14, true);
    test_sites2(sites, 36, "7mer", 36, 15, true);
    test_sites2(sites, 37, "7mer", 37, 18, true);
    test_sites2(sites, 38, "7mer", 38, 19, true);
    test_sites2(sites, 39, "7mer", 39, 20, true);
    test_sites2(sites, 40, "7mer", 40, 21, true);
    test_sites2(sites, 41, "7mer", 41, 32, true);
    test_sites2(sites, 42, "7mer", 42, 33, true);
    test_sites2(sites, 43, "7mer", 43, 34, true);

    test_sites2(sites, 44, "7mer", 44, 2, true);
    test_sites2(sites, 45, "7mer", 45, 13, true);
    test_sites2(sites, 46, "7mer", 46, 14, true);
    test_sites2(sites, 47, "7mer", 47, 15, true);
    test_sites2(sites, 48, "7mer", 48, 18, true);
    test_sites2(sites, 49, "7mer", 49, 19, true);
    test_sites2(sites, 50, "7mer", 50, 20, true);
    test_sites2(sites, 51, "7mer", 51, 21, true);
    test_sites2(sites, 52, "7mer", 52, 32, true);
    test_sites2(sites, 53, "7mer", 53, 33, true);
    test_sites2(sites, 54, "7mer", 54, 34, true);
}

TEST_F(Site01Nmer1, mir124_def) {
    read_files(false);
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef1 = "7mGU+";
    int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef1, mOverlapDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(55u, sites.get_length());

    test_sites2(sites, 0, "", 0, 0, false);
    test_sites2(sites, 1, "7mer", 1, 1, false);
    test_sites2(sites, 2, "7mer", 2, 13, false);
    test_sites2(sites, 3, "7mer", 3, 14, false);
    test_sites2(sites, 4, "7mer", 4, 15, false);
    test_sites2(sites, 5, "7mer", 5, 18, false);
    test_sites2(sites, 6, "7mer", 6, 19, false);
    test_sites2(sites, 7, "7mer", 7, 20, false);
    test_sites2(sites, 8, "7mer", 8, 32, false);
    test_sites2(sites, 9, "7mer", 9, 33, false);
    test_sites2(sites, 10, "7mer", 10, 34, false);

    test_sites2(sites, 11, "", 11, 0, false);
    test_sites2(sites, 12, "7mer", 12, 1, false);
    test_sites2(sites, 13, "7mer", 13, 13, false);
    test_sites2(sites, 14, "7mer", 14, 14, false);
    test_sites2(sites, 15, "7mer", 15, 15, false);
    test_sites2(sites, 16, "7mer", 16, 18, false);
    test_sites2(sites, 17, "7mer", 17, 19, false);
    test_sites2(sites, 18, "7mer", 18, 20, false);
    test_sites2(sites, 19, "7mer", 19, 32, false);
    test_sites2(sites, 20, "7mer", 20, 33, false);
    test_sites2(sites, 21, "7mer", 21, 34, false);

    test_sites2(sites, 22, "7mer", 22, 1, true);
    test_sites2(sites, 23, "7mer", 23, 13, true);
    test_sites2(sites, 24, "7mer", 24, 14, true);
    test_sites2(sites, 25, "7mer", 25, 15, true);
    test_sites2(sites, 26, "7mer", 26, 18, true);
    test_sites2(sites, 27, "7mer", 27, 19, true);
    test_sites2(sites, 28, "7mer", 28, 20, true);
    test_sites2(sites, 29, "7mer", 29, 21, true);
    test_sites2(sites, 30, "7mer", 30, 32, true);
    test_sites2(sites, 31, "7mer", 31, 33, true);
    test_sites2(sites, 32, "7mer", 32, 34, true);

    test_sites2(sites, 33, "7mer", 33, 1, true);
    test_sites2(sites, 34, "7mer", 34, 13, true);
    test_sites2(sites, 35, "7mer", 35, 14, true);
    test_sites2(sites, 36, "7mer", 36, 15, true);
    test_sites2(sites, 37, "7mer", 37, 18, true);
    test_sites2(sites, 38, "7mer", 38, 19, true);
    test_sites2(sites, 39, "7mer", 39, 20, true);
    test_sites2(sites, 40, "7mer", 40, 21, true);
    test_sites2(sites, 41, "7mer", 41, 32, true);
    test_sites2(sites, 42, "7mer", 42, 33, true);
    test_sites2(sites, 43, "7mer", 43, 34, true);

    test_sites2(sites, 44, "7mer", 44, 2, true);
    test_sites2(sites, 45, "7mer", 45, 13, true);
    test_sites2(sites, 46, "7mer", 46, 14, true);
    test_sites2(sites, 47, "7mer", 47, 15, true);
    test_sites2(sites, 48, "7mer", 48, 18, true);
    test_sites2(sites, 49, "7mer", 49, 19, true);
    test_sites2(sites, 50, "7mer", 50, 20, true);
    test_sites2(sites, 51, "7mer", 51, 21, true);
    test_sites2(sites, 52, "7mer", 52, 32, true);
    test_sites2(sites, 53, "7mer", 53, 33, true);
    test_sites2(sites, 54, "7mer", 54, 34, true);
}

}
