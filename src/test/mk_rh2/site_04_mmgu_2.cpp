#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site04MMGU2 : public TestSiteRH2 {
protected:
    Site04MMGU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mGU1";
        mOverlapDef = "orig";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef rh2mfe::RH2SeedSites TSit;
    typedef rh2mfe::RH2SeedSeqs TSeed;

};

TEST_F(Site04MMGU2, mir1_mmgu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    seedSeqs.set_mirna_seq(mirna_seqs[1]);
    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs();

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(26u, sites.get_length());

    test_sites2(sites, 0, "6mer", 0, 24, true);
    test_sites2(sites, 1, "6mer", 1, 24, true);

    test_sites2(sites, 2, "6mer_GUT", 2, 24, true);
    test_sites2(sites, 3, "6mer_GUT", 3, 24, true);
    test_sites2(sites, 4, "6mer_GUT", 4, 24, true);
    test_sites2(sites, 5, "6mer_GUT", 5, 24, true);
    test_sites2(sites, 6, "6mer_GUT", 6, 24, true);
    test_sites2(sites, 7, "6mer_GUT", 7, 24, true);

    test_sites2(sites, 8, "6mer_GUM", 8, 24, true);
    test_sites2(sites, 9, "6mer_GUM", 9, 24, true);
    test_sites2(sites, 10, "6mer_GUM", 10, 24, true);
    test_sites2(sites, 11, "6mer_GUM", 11, 24, true);
    test_sites2(sites, 12, "6mer_GUM", 12, 24, true);
    test_sites2(sites, 13, "6mer_GUM", 13, 24, true);

    test_sites2(sites, 14, "6mer_GUT", 14, 24, true);
    test_sites2(sites, 15, "6mer_GUT", 15, 24, true);
    test_sites2(sites, 16, "6mer_GUT", 16, 24, true);
    test_sites2(sites, 17, "6mer_GUT", 17, 24, true);
    test_sites2(sites, 18, "6mer_GUT", 18, 24, true);
    test_sites2(sites, 19, "6mer_GUT", 19, 24, true);

    test_sites2(sites, 20, "6mer_GUT", 20, 24, true);
    test_sites2(sites, 21, "6mer_GUT", 21, 24, true);
    test_sites2(sites, 22, "6mer_GUT", 22, 24, true);
    test_sites2(sites, 23, "6mer_GUT", 23, 24, true);
    test_sites2(sites, 24, "6mer_GUT", 24, 24, true);
    test_sites2(sites, 25, "6mer_GUT", 25, 24, true);
}

TEST_F(Site04MMGU2, mir1_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    mSeedDef[0] = "7mGU+";
    seedSeqs.set_mirna_seq(mirna_seqs[1]);
    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs();

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(104u, sites.get_length());

    test_sites2(sites, 0, "7mer_GUM", 0, 24, true);
    test_sites2(sites, 1, "7mer_GUM", 1, 24, true);
//    test_sites2(sites, 2, "7mer_GUT", 2, 24, false);
//    test_sites2(sites, 3, "7mer_GUT", 3, 24, false);
    test_sites2(sites, 2, "7mer_GUT", 4, 24, true);
    test_sites2(sites, 3, "7mer_GUT", 5, 24, true);
//    test_sites2(sites, 6, "7mer_GUT", 6, 24, false);
//    test_sites2(sites, 7, "7mer_GUT", 7, 24, false);
//    test_sites2(sites, 8, "7mer_GUM", 8, 24, false);
//    test_sites2(sites, 9, "7mer_GUM", 9, 24, false);
    test_sites2(sites, 4, "7mer_GUM", 10, 24, true);
    test_sites2(sites, 5, "7mer_GUM", 11, 24, true);
//    test_sites2(sites, 12, "7mer_GUM", 12, 24, false);
//    test_sites2(sites, 13, "7mer_GUM", 13, 24, false);
//    test_sites2(sites, 14, "7mer_GUT", 14, 24, false);
//    test_sites2(sites, 15, "7mer_GUT", 15, 24, false);
    test_sites2(sites, 6, "7mer_GUT", 16, 24, true);
    test_sites2(sites, 7, "7mer_GUT", 17, 24, true);
//    test_sites2(sites, 18, "7mer_GUT", 18, 24, false);
//    test_sites2(sites, 19, "7mer_GUT", 19, 24, false);
//    test_sites2(sites, 20, "7mer_GUT", 20, 24, false);
//    test_sites2(sites, 21, "7mer_GUT", 21, 24, false);
    test_sites2(sites, 8, "7mer_GUT", 22, 24, true);
    test_sites2(sites, 9, "7mer_GUT", 23, 24, true);
//    test_sites2(sites, 24, "7mer_GUT", 24, 24, false);
//    test_sites2(sites, 25, "7mer_GUT", 25, 24, false);
    test_sites2(sites, 10, "7mer_GUT", 0, 14, true);
    test_sites2(sites, 11, "7mer_GUT", 1, 14, true);
    test_sites2(sites, 12, "7mer_GUT", 2, 14, true);
    test_sites2(sites, 13, "7mer_GUT", 3, 14, true);
    test_sites2(sites, 14, "7mer_GUT", 4, 14, true);
    test_sites2(sites, 15, "7mer_GUT", 5, 14, true);
    test_sites2(sites, 16, "7mer_GUT", 6, 14, true);
    test_sites2(sites, 17, "7mer_GUT", 7, 14, true);
    test_sites2(sites, 18, "7mer_GUT", 8, 14, true);
    test_sites2(sites, 19, "7mer_GUT", 9, 14, true);
    test_sites2(sites, 20, "7mer_GUT", 10, 14, true);
    test_sites2(sites, 21, "7mer_GUT", 11, 14, true);
    test_sites2(sites, 22, "7mer_GUT", 12, 14, true);
    test_sites2(sites, 23, "7mer_GUT", 13, 14, true);
    test_sites2(sites, 24, "7mer_GUT", 14, 14, true);
    test_sites2(sites, 25, "7mer_GUT", 15, 14, true);
    test_sites2(sites, 26, "7mer_GUT", 16, 14, true);
    test_sites2(sites, 27, "7mer_GUT", 17, 14, true);
    test_sites2(sites, 28, "7mer_GUT", 18, 14, true);
    test_sites2(sites, 29, "7mer_GUT", 19, 14, true);
    test_sites2(sites, 30, "7mer_GUT", 20, 14, true);
    test_sites2(sites, 31, "7mer_GUT", 21, 14, true);
    test_sites2(sites, 32, "7mer_GUT", 22, 14, true);
    test_sites2(sites, 33, "7mer_GUT", 23, 14, true);
    test_sites2(sites, 34, "7mer_GUT", 24, 14, true);
    test_sites2(sites, 35, "7mer_GUT", 25, 14, true);
    test_sites2(sites, 36, "7mer_GUT", 26, 14, true);
    test_sites2(sites, 37, "7mer_GUT", 27, 14, true);
    test_sites2(sites, 38, "7mer_GUT", 28, 14, true);
    test_sites2(sites, 39, "7mer_GUT", 29, 14, true);
    test_sites2(sites, 40, "7mer_GUT", 30, 14, true);
    test_sites2(sites, 41, "7mer_GUT", 31, 14, true);
    test_sites2(sites, 42, "7mer_GUT", 32, 14, true);
    test_sites2(sites, 43, "7mer_GUT", 33, 14, true);
    test_sites2(sites, 44, "7mer_GUT", 34, 14, true);
    test_sites2(sites, 45, "7mer_GUT", 35, 14, true);
    test_sites2(sites, 46, "7mer_GUT", 36, 14, true);
    test_sites2(sites, 47, "7mer_GUT", 37, 14, true);
    test_sites2(sites, 48, "7mer_GUT", 38, 14, true);
    test_sites2(sites, 49, "7mer_GUT", 39, 14, true);
    test_sites2(sites, 50, "7mer_GUT", 40, 14, true);
    test_sites2(sites, 51, "7mer_GUT", 41, 14, true);
    test_sites2(sites, 52, "7mer_GUT", 42, 14, true);
    test_sites2(sites, 53, "7mer_GUT", 43, 14, true);
    test_sites2(sites, 54, "7mer_GUT", 44, 14, true);
    test_sites2(sites, 55, "7mer_GUT", 45, 14, true);
    test_sites2(sites, 56, "7mer_GUT", 46, 14, true);
    test_sites2(sites, 57, "7mer_GUT", 47, 14, true);
    test_sites2(sites, 58, "7mer_GUT", 48, 14, true);
    test_sites2(sites, 59, "7mer_GUT", 49, 14, true);
    test_sites2(sites, 60, "7mer_GUT", 50, 14, true);
    test_sites2(sites, 61, "7mer_GUT", 51, 14, true);
    test_sites2(sites, 62, "7mer_GUT", 52, 14, true);
    test_sites2(sites, 63, "7mer_GUT", 53, 14, true);
    test_sites2(sites, 64, "7mer_GUT", 54, 14, true);
    test_sites2(sites, 65, "7mer_GUT", 55, 14, true);
    test_sites2(sites, 66, "7mer_GUT", 56, 14, true);
    test_sites2(sites, 67, "7mer_GUT", 57, 14, true);
    test_sites2(sites, 68, "7mer_GUT", 58, 14, true);
    test_sites2(sites, 69, "7mer_GUT", 59, 14, true);
    test_sites2(sites, 70, "7mer_GUT", 60, 14, true);
    test_sites2(sites, 71, "7mer_GUT", 61, 14, true);
    test_sites2(sites, 72, "7mer_GUT", 62, 14, true);
    test_sites2(sites, 73, "7mer_GUT", 63, 14, true);
    test_sites2(sites, 74, "7mer_GUT", 64, 14, true);
    test_sites2(sites, 75, "7mer_GUT", 65, 14, true);
    test_sites2(sites, 76, "7mer_GUT", 66, 14, true);
    test_sites2(sites, 77, "7mer_GUT", 67, 14, true);
    test_sites2(sites, 78, "7mer_GUT", 68, 14, true);
    test_sites2(sites, 79, "7mer_GUT", 69, 14, true);
    test_sites2(sites, 80, "7mer_GUT", 70, 14, true);
    test_sites2(sites, 81, "7mer_GUT", 71, 14, true);
    test_sites2(sites, 82, "7mer_GUT", 72, 14, true);
    test_sites2(sites, 83, "7mer_GUT", 73, 14, true);
    test_sites2(sites, 84, "7mer_GUT", 74, 14, true);
    test_sites2(sites, 85, "7mer_GUT", 75, 14, true);
    test_sites2(sites, 86, "7mer_GUT", 76, 14, true);
    test_sites2(sites, 87, "7mer_GUT", 77, 14, true);
    test_sites2(sites, 88, "7mer_GUT", 78, 14, true);
    test_sites2(sites, 89, "7mer_GUT", 79, 14, true);
    test_sites2(sites, 90, "7mer_GUT", 80, 14, true);
    test_sites2(sites, 91, "7mer_GUT", 81, 14, true);
    test_sites2(sites, 92, "7mer_GUT", 82, 14, true);
    test_sites2(sites, 93, "7mer_GUT", 83, 14, true);
    test_sites2(sites, 94, "7mer_GUT", 84, 14, true);
    test_sites2(sites, 95, "7mer_GUT", 85, 14, true);
    test_sites2(sites, 96, "7mer_GUT", 86, 14, true);
    test_sites2(sites, 97, "7mer_GUT", 87, 14, true);
    test_sites2(sites, 98, "7mer_GUT", 88, 14, true);
    test_sites2(sites, 99, "7mer_GUT", 89, 14, true);
    test_sites2(sites, 100, "7mer_GUT", 90, 14, true);
    test_sites2(sites, 101, "7mer_GUT", 91, 14, true);
    test_sites2(sites, 102, "7mer_GUT", 92, 14, true);
    test_sites2(sites, 103, "7mer_GUT", 93, 14, true);
}

}
