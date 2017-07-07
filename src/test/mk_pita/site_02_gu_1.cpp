#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site02GU1 : public TestSitePITA {
protected:
    Site02GU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_pita/";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "1";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "0";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef ptddg::PITASeedSites<mikan::TRNATYPE> TSit;

};

TEST_F(Site02GU1, mir124_gu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(11u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);

    test_sites(sites, 1, "GUT", 1, 24, false, 0);
    test_sites(sites, 2, "7mer_GUT", 2, 24, true, 2);
    test_sites(sites, 3, "7mer_GUT", 3, 24, true, 2);
    test_sites(sites, 4, "8mer_GUT", 4, 24, true, 2);
    test_sites(sites, 5, "8mer_GUT", 5, 24, true, 2);

    test_sites(sites, 6, "GUT", 6, 24, false, 0);
    test_sites(sites, 7, "7mer_GUT", 7, 24, true, 3);
    test_sites(sites, 8, "7mer_GUT", 8, 24, true, 3);
    test_sites(sites, 9, "8mer_GUT", 9, 24, true, 3);
    test_sites(sites, 10, "8mer_GUT", 10, 24, true, 3);
}

TEST_F(Site02GU1, mir124_gu_plus) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "+";
    int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(16u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);

    test_sites(sites, 1, "GU+", 1, 24, false, 0);
    test_sites(sites, 2, "7mer_GU+", 2, 24, true, 0);
    test_sites(sites, 3, "7mer_GU+", 3, 24, true, 0);
    test_sites(sites, 4, "8mer_GU+", 4, 24, true, 0);
    test_sites(sites, 5, "8mer_GU+", 5, 24, true, 0);

    test_sites(sites, 6, "GU+", 6, 24, false, 0);
    test_sites(sites, 7, "7mer_GU+", 7, 24, true, 0);
    test_sites(sites, 8, "7mer_GU+", 8, 24, true, 0);
    test_sites(sites, 9, "8mer_GU+", 9, 24, true, 0);
    test_sites(sites, 10, "8mer_GU+", 10, 24, true, 0);

    test_sites(sites, 11, "GU+", 11, 24, false, 0);
    test_sites(sites, 12, "7mer_GU+", 12, 24, true, 0);
    test_sites(sites, 13, "7mer_GU+", 13, 24, true, 0);
    test_sites(sites, 14, "8mer_GU+", 14, 24, true, 0);
    test_sites(sites, 15, "8mer_GU+", 15, 24, true, 0);
}

TEST_F(Site02GU1, mir124_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    int ret_val = sites.find_seed_sites(mirna_seqs[0], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(11u, sites.get_length());

    test_sites(sites, 0, "7mer", 0, 24, true, 0);

    test_sites(sites, 1, "GUT", 1, 24, false, 0);
    test_sites(sites, 2, "7mer_GUT", 2, 24, true, 2);
    test_sites(sites, 3, "7mer_GUT", 3, 24, true, 2);
    test_sites(sites, 4, "8mer_GUT", 4, 24, true, 2);
    test_sites(sites, 5, "8mer_GUT", 5, 24, true, 2);

    test_sites(sites, 6, "GUT", 6, 24, false, 0);
    test_sites(sites, 7, "7mer_GUT", 7, 24, true, 3);
    test_sites(sites, 8, "7mer_GUT", 8, 24, true, 3);
    test_sites(sites, 9, "8mer_GUT", 9, 24, true, 3);
    test_sites(sites, 10, "8mer_GUT", 10, 24, true, 3);

}
}
