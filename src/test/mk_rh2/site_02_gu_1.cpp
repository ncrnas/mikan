#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site02GU1 : public TestSiteRH2 {
protected:
    Site02GU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_1.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mGU1";
        mOverlapDef = "orig";
    }

};

TEST_F(Site02GU1, mir124_6mer_gu1) {
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
    EXPECT_EQ(11u, sites.get_length());

    test_sites2(sites, 0, "6mer", 0, 24, true);

    test_sites2(sites, 1, "6mer_GUT", 1, 24, true);
    test_sites2(sites, 2, "6mer_GUT", 2, 24, true);
    test_sites2(sites, 3, "6mer_GUT", 3, 24, true);
    test_sites2(sites, 4, "6mer_GUT", 4, 24, true);
    test_sites2(sites, 5, "6mer_GUT", 5, 24, true);

    test_sites2(sites, 6, "6mer_GUT", 6, 24, true);
    test_sites2(sites, 7, "6mer_GUT", 7, 24, true);
    test_sites2(sites, 8, "6mer_GUT", 8, 24, true);
    test_sites2(sites, 9, "6mer_GUT", 9, 24, true);
    test_sites2(sites, 10, "6mer_GUT", 10, 24, true);
}

TEST_F(Site02GU1, mir124_6mer_gu_plus) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TOp ops;
    TSeed seedSeqs(ops);

    mSeedDef[0] = "6mGU+";

    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(16u, sites.get_length());

    test_sites2(sites, 0, "6mer", 0, 24, true);

    test_sites2(sites, 1, "6mer_GUT", 1, 24, true);
    test_sites2(sites, 2, "6mer_GUT", 2, 24, true);
    test_sites2(sites, 3, "6mer_GUT", 3, 24, true);
    test_sites2(sites, 4, "6mer_GUT", 4, 24, true);
    test_sites2(sites, 5, "6mer_GUT", 5, 24, true);

    test_sites2(sites, 6, "6mer_GUT", 6, 24, true);
    test_sites2(sites, 7, "6mer_GUT", 7, 24, true);
    test_sites2(sites, 8, "6mer_GUT", 8, 24, true);
    test_sites2(sites, 9, "6mer_GUT", 9, 24, true);
    test_sites2(sites, 10, "6mer_GUT", 10, 24, true);

    test_sites2(sites, 11, "6mer_GUT", 11, 24, true);
    test_sites2(sites, 12, "6mer_GUT", 12, 24, true);
    test_sites2(sites, 13, "6mer_GUT", 13, 24, true);
    test_sites2(sites, 14, "6mer_GUT", 14, 24, true);
    test_sites2(sites, 15, "6mer_GUT", 15, 24, true);
}

TEST_F(Site02GU1, mir124_7mer_gu1) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TOp ops;
    TSeed seedSeqs(ops);

    mSeedDef[0] = "7mGU1";

    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(9u, sites.get_length());

    test_sites2(sites, 0, "7mer", 0, 24, true);
//    test_sites2(sites, 1, "7mer_GUT", 1, 24, false);
    test_sites2(sites, 1, "7mer_GUT", 2, 24, true);
    test_sites2(sites, 2, "7mer_GUT", 3, 24, true);
    test_sites2(sites, 3, "7mer_GUT", 4, 24, true);
    test_sites2(sites, 4, "7mer_GUT", 5, 24, true);
//    test_sites2(sites, 6, "7mer_GUT", 6, 24, false);
    test_sites2(sites, 5, "7mer_GUT", 7, 24, true);
    test_sites2(sites, 6, "7mer_GUT", 8, 24, true);
    test_sites2(sites, 7, "7mer_GUT", 9, 24, true);
    test_sites2(sites, 8, "7mer_GUT", 10, 24, true);
}

TEST_F(Site02GU1, mir124_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TOp ops;
    TSeed seedSeqs(ops);

    mSeedDef[0] = "7mGU+";

    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(13u, sites.get_length());

    test_sites2(sites, 0, "7mer", 0, 24, true);
//    test_sites2(sites, 1, "7mer_GUT", 1, 24, false);
    test_sites2(sites, 1, "7mer_GUT", 2, 24, true);
    test_sites2(sites, 2, "7mer_GUT", 3, 24, true);
    test_sites2(sites, 3, "7mer_GUT", 4, 24, true);
    test_sites2(sites, 4, "7mer_GUT", 5, 24, true);
//    test_sites2(sites, 6, "7mer_GUT", 6, 24, false);
    test_sites2(sites, 5, "7mer_GUT", 7, 24, true);
    test_sites2(sites, 6, "7mer_GUT", 8, 24, true);
    test_sites2(sites, 7, "7mer_GUT", 9, 24, true);
    test_sites2(sites, 8, "7mer_GUT", 10, 24, true);
//    test_sites2(sites, 11, "7mer_GUT", 11, 24, false);
    test_sites2(sites, 9, "7mer_GUT", 12, 24, true);
    test_sites2(sites, 10, "7mer_GUT", 13, 24, true);
    test_sites2(sites, 11, "7mer_GUT", 14, 24, true);
    test_sites2(sites, 12, "7mer_GUT", 15, 24, true);
}
}
