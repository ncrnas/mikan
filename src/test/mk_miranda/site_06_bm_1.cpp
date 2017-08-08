#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_miranda.hpp"

namespace {

class Site06BM1 : public TestSiteMR3AS {
protected:
    Site06BM1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_miranda/";

        resize(mSeedDef, 6);
        mSeedDef[0] = 'Y';
        mSeedDef[1] = 'Y';
        mSeedDef[2] = 'Y';
        mSeedDef[3] = "1";
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "1";
    }

};

TEST_F(Site06BM1, mir124_bm) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    TOp ops;
    TSeed seedSeqs(ops);

    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(4u, sites.get_length());

//    test_sites(sites, 0, "MM", 0, 25, false, 0);
//    test_sites(sites, 1, "MM", 1, 25, false, 0);
    test_sites(sites, 0, "8mer_MM", 6, 26, true, 5);
    test_sites(sites, 1, "8mer_MM", 7, 26, true, 5);
    test_sites(sites, 2, "8mer_MMGU", 4, 26, true, 5);
    test_sites(sites, 3, "8mer_MMGU", 5, 26, true, 5);

//    test_sites(sites, 6, "BT", 6, 25, false, 0);
//    test_sites(sites, 7, "BT", 7, 25, false, 0);
//    test_sites(sites, 8, "BT", 7, 26, false, 0);
//    test_sites(sites, 9, "BT", 6, 26, false, 0);
}

TEST_F(Site06BM1, mir124_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "+";
    mSeedDef[4] = "1:1";
    mSeedDef[5] = "1";
    TOp ops;
    TSeed seedSeqs(ops);

    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(4u, sites.get_length());

//    test_sites(sites, 0, "MM", 0, 25, false, 0);
//    test_sites(sites, 1, "MM", 1, 25, false, 0);
    test_sites(sites, 0, "8mer_MM", 6, 26, true, 5);
    test_sites(sites, 1, "8mer_MM", 7, 26, true, 5);
    test_sites(sites, 2, "8mer_MMGU", 4, 26, true, 5);
    test_sites(sites, 3, "8mer_MMGU", 5, 26, true, 5);

//    test_sites(sites, 6, "BT", 6, 25, false, 0);
//    test_sites(sites, 7, "BT", 7, 25, false, 0);
//    test_sites(sites, 8, "BT", 7, 26, false, 0);
//    test_sites(sites, 9, "BT", 6, 26, false, 0);
}
}
