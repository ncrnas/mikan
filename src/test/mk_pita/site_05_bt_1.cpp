#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site05BT1 : public TestSitePITA {
protected:
    Site05BT1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_pita/";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "0";
        mSeedDef[4] = "0:0";
        mSeedDef[5] = "1";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef ptddg::PITASeedSites TSit;
    typedef ptddg::PITASeedSeqs TSeed;

};

TEST_F(Site05BT1, mir124_bt) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    seedSeqs.set_mirna_seq(mirna_seqs[0]);
    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs();

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(0u, sites.get_length());
}

TEST_F(Site05BT1, mir124_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    seedSeqs.set_mirna_seq(mirna_seqs[0]);
    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs();

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(18u, sites.get_length());

    test_sites(sites, 0, "MM", 4, 25, false, 0);
    test_sites(sites, 1, "MM", 5, 25, false, 0);

    test_sites(sites, 2, "MM", 0, 25, false, 0);
    test_sites(sites, 3, "MM", 1, 25, false, 0);
    test_sites(sites, 4, "MM", 2, 25, false, 0);
    test_sites(sites, 5, "MM", 3, 25, false, 0);

    test_sites(sites, 6, "8mer_MM", 9, 24, true, 4);

    test_sites(sites, 7, "8mer_MM", 19, 24, true, 4);

    test_sites(sites, 8, "MM", 16, 24, false, 0);

    test_sites(sites, 9, "8mer_MM", 17, 24, true, 4);
    test_sites(sites, 10, "8mer_MM", 18, 24, true, 4);

    test_sites(sites, 11, "8mer_MM", 23, 24, true, 5);
    test_sites(sites, 12, "8mer_MM", 24, 24, true, 5);

    test_sites(sites, 13, "MM", 20, 24, false, 0);
    test_sites(sites, 14, "8mer_MM", 21, 24, true, 5);
    test_sites(sites, 15, "8mer_MM", 22, 24, true, 5);

    test_sites(sites, 16, "8mer_MMGU", 10, 24, true, 4);

    test_sites(sites, 17, "8mer_MMGU", 15, 24, true, 4);
}

}
