#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site01Nmer1 : public TestSitePITA {
protected:
    Site01Nmer1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";
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
        mSeedDef[5] = "0";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef ptddg::PITASeedSites TSit;
    typedef ptddg::PITASeedSeqs TSeed;

};

TEST_F(Site01Nmer1, mir124_8mer) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(28u, sites.get_length());

//    test_sites(sites, 0, "", 0, 0, false, 0);
//    test_sites(sites, 1, "", 1, 1, false, 0);
//    test_sites(sites, 2, "", 2, 13, false, 0);
//    test_sites(sites, 3, "", 3, 14, false, 0);
//    test_sites(sites, 4, "", 4, 15, false, 0);
    test_sites(sites, 0, "6mer", 5, 18, true, 0);
    test_sites(sites, 1, "6mer", 6, 19, true, 0);
    test_sites(sites, 2, "6mer", 7, 20, true, 0);
    test_sites(sites, 3, "6mer", 8, 32, true, 0);
    test_sites(sites, 4, "6mer", 9, 33, true, 0);
//    test_sites(sites, 10, "", 10, 34, false, 0);

//    test_sites(sites, 11, "", 11, 0, false, 0);
//    test_sites(sites, 12, "", 12, 1, false, 0);
//    test_sites(sites, 13, "", 13, 13, false, 0);
//    test_sites(sites, 14, "", 14, 14, false, 0);
//    test_sites(sites, 15, "", 15, 15, false, 0);
    test_sites(sites, 5, "6mer", 16, 18, true, 0);
    test_sites(sites, 6, "6mer", 17, 19, true, 0);
    test_sites(sites, 7, "6mer", 18, 20, true, 0);
    test_sites(sites, 8, "6mer", 19, 32, true, 0);
    test_sites(sites, 9, "6mer", 20, 33, true, 0);
//    test_sites(sites, 21, "", 21, 34, false, 0);

//    test_sites(sites, 22, "", 22, 1, false, 0);
//    test_sites(sites, 23, "", 23, 13, false, 0);
//    test_sites(sites, 24, "", 24, 14, false, 0);
//    test_sites(sites, 25, "", 25, 15, false, 0);
    test_sites(sites, 10, "7mer", 26, 18, true, 0);
    test_sites(sites, 11, "7mer", 27, 19, true, 0);
    test_sites(sites, 12, "7mer", 28, 20, true, 0);
    test_sites(sites, 13, "7mer", 29, 21, true, 0);
    test_sites(sites, 14, "7mer", 30, 32, true, 0);
    test_sites(sites, 15, "7mer", 31, 33, true, 0);
//    test_sites(sites, 32, "", 32, 34, false, 0);

//    test_sites(sites, 33, "", 33, 1, false, 0);
//    test_sites(sites, 34, "", 34, 13, false, 0);
//    test_sites(sites, 35, "", 35, 14, false, 0);
//    test_sites(sites, 36, "", 36, 15, false, 0);
    test_sites(sites, 16, "7mer", 37, 18, true, 0);
    test_sites(sites, 17, "7mer", 38, 19, true, 0);
    test_sites(sites, 18, "7mer", 39, 20, true, 0);
    test_sites(sites, 19, "7mer", 40, 21, true, 0);
    test_sites(sites, 20, "7mer", 41, 32, true, 0);
    test_sites(sites, 21, "7mer", 42, 33, true, 0);
//    test_sites(sites, 43, "", 43, 34, false, 0);

//    test_sites(sites, 44, "", 44, 2, false, 0);
//    test_sites(sites, 45, "", 45, 13, false, 0);
//    test_sites(sites, 46, "", 46, 14, false, 0);
//    test_sites(sites, 47, "", 47, 15, false, 0);
    test_sites(sites, 22, "8mer", 48, 18, true, 0);
    test_sites(sites, 23, "8mer", 49, 19, true, 0);
    test_sites(sites, 24, "8mer", 50, 20, true, 0);
    test_sites(sites, 25, "8mer", 51, 21, true, 0);
    test_sites(sites, 26, "8mer", 52, 32, true, 0);
    test_sites(sites, 27, "8mer", 53, 33, true, 0);
//    test_sites(sites, 54, "", 54, 34, false, 0);
}

TEST_F(Site01Nmer1, mir124_7mer) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    mSeedDef[2] = 'N';

    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(28u, sites.get_length());

    test_sites(sites, 3, "6mer", 8, 32, true, 0);
    test_sites(sites, 8, "6mer", 19, 32, true, 0);
    test_sites(sites, 14, "7mer", 30, 32, true, 0);
    test_sites(sites, 20, "7mer", 41, 32, true, 0);
    test_sites(sites, 26, "7mer", 52, 32, true, 0);
}

TEST_F(Site01Nmer1, mir124_6mer) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    mSeedDef[1] = 'N';
    mSeedDef[2] = 'N';

    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(28u, sites.get_length());

    test_sites(sites, 3, "6mer", 8, 32, true, 0);
    test_sites(sites, 8, "6mer", 19, 32, true, 0);
    test_sites(sites, 14, "6mer", 30, 32, true, 0);
    test_sites(sites, 20, "6mer", 41, 32, true, 0);
    test_sites(sites, 26, "6mer", 52, 32, true, 0);
}

TEST_F(Site01Nmer1, mir124_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";

    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs(mirna_seqs[0]);

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(28u, sites.get_length());

//    test_sites(sites, 0, "", 0, 0, false, 0);
//    test_sites(sites, 1, "", 1, 1, false, 0);
//    test_sites(sites, 2, "", 2, 13, false, 0);
//    test_sites(sites, 3, "", 3, 14, false, 0);
//    test_sites(sites, 4, "", 4, 15, false, 0);
    test_sites(sites, 0, "6mer", 5, 18, true, 0);
    test_sites(sites, 1, "6mer", 6, 19, true, 0);
    test_sites(sites, 2, "6mer", 7, 20, true, 0);
    test_sites(sites, 3, "6mer", 8, 32, true, 0);
    test_sites(sites, 4, "6mer", 9, 33, true, 0);
//    test_sites(sites, 10, "", 10, 34, false, 0);

//    test_sites(sites, 11, "", 11, 0, false, 0);
//    test_sites(sites, 12, "", 12, 1, false, 0);
//    test_sites(sites, 13, "", 13, 13, false, 0);
//    test_sites(sites, 14, "", 14, 14, false, 0);
//    test_sites(sites, 15, "", 15, 15, false, 0);
    test_sites(sites, 5, "6mer", 16, 18, true, 0);
    test_sites(sites, 6, "6mer", 17, 19, true, 0);
    test_sites(sites, 7, "6mer", 18, 20, true, 0);
    test_sites(sites, 8, "6mer", 19, 32, true, 0);
    test_sites(sites, 9, "6mer", 20, 33, true, 0);
//    test_sites(sites, 21, "", 21, 34, false, 0);

//    test_sites(sites, 22, "", 22, 1, false, 0);
//    test_sites(sites, 23, "", 23, 13, false, 0);
//    test_sites(sites, 24, "", 24, 14, false, 0);
//    test_sites(sites, 25, "", 25, 15, false, 0);
    test_sites(sites, 10, "7mer", 26, 18, true, 0);
    test_sites(sites, 11, "7mer", 27, 19, true, 0);
    test_sites(sites, 12, "7mer", 28, 20, true, 0);
    test_sites(sites, 13, "7mer", 29, 21, true, 0);
    test_sites(sites, 14, "7mer", 30, 32, true, 0);
    test_sites(sites, 15, "7mer", 31, 33, true, 0);
//    test_sites(sites, 32, "", 32, 34, false, 0);

//    test_sites(sites, 33, "", 33, 1, false, 0);
//    test_sites(sites, 34, "", 34, 13, false, 0);
//    test_sites(sites, 35, "", 35, 14, false, 0);
//    test_sites(sites, 36, "", 36, 15, false, 0);
    test_sites(sites, 16, "7mer", 37, 18, true, 0);
    test_sites(sites, 17, "7mer", 38, 19, true, 0);
    test_sites(sites, 18, "7mer", 39, 20, true, 0);
    test_sites(sites, 19, "7mer", 40, 21, true, 0);
    test_sites(sites, 20, "7mer", 41, 32, true, 0);
    test_sites(sites, 21, "7mer", 42, 33, true, 0);
//    test_sites(sites, 43, "", 43, 34, false, 0);

//    test_sites(sites, 44, "", 44, 2, false, 0);
//    test_sites(sites, 45, "", 45, 13, false, 0);
//    test_sites(sites, 46, "", 46, 14, false, 0);
//    test_sites(sites, 47, "", 47, 15, false, 0);
    test_sites(sites, 22, "8mer", 48, 18, true, 0);
    test_sites(sites, 23, "8mer", 49, 19, true, 0);
    test_sites(sites, 24, "8mer", 50, 20, true, 0);
    test_sites(sites, 25, "8mer", 51, 21, true, 0);
    test_sites(sites, 26, "8mer", 52, 32, true, 0);
    test_sites(sites, 27, "8mer", 53, 33, true, 0);
//    test_sites(sites, 54, "", 54, 34, false, 0);
}

}
