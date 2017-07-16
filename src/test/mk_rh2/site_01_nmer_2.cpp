#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site01Nmer2 : public TestSiteRH2 {
protected:
    Site01Nmer2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        resize(mSeedDef, 1);
        mSeedDef[0] = "6mer";
        mOverlapDef = "orig";
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef rh2mfe::RH2SeedSites TSit;
    typedef rh2mfe::RH2SeedSeqs TSeed;

};

TEST_F(Site01Nmer2, mir1_6mer) {
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
    EXPECT_EQ(5u, sites.get_length());

    test_sites2(sites, 0, "6mer", 0, 25, true);
    test_sites2(sites, 1, "6mer", 1, 25, true);
    test_sites2(sites, 2, "6mer", 2, 25, true);
    test_sites2(sites, 3, "6mer", 3, 25, true);
    test_sites2(sites, 4, "6mer", 4, 25, true);
}

TEST_F(Site01Nmer2, mir1_7mer) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    TSeed seedSeqs;

    mSeedDef[0] = "7mer";
    seedSeqs.set_mirna_seq(mirna_seqs[1]);
    seedSeqs.set_flags(mSeedDef);
    seedSeqs.create_seed_seqs();

    int ret_val = sites.find_seed_sites(seedSeqs, mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(5u, sites.get_length());

    test_sites2(sites, 0, "7mer", 0, 25, false);
    test_sites2(sites, 1, "7mer", 1, 25, false);
    test_sites2(sites, 2, "7mer", 2, 25, true);
    test_sites2(sites, 3, "7mer", 3, 25, true);
    test_sites2(sites, 4, "7mer", 4, 25, true);
}

TEST_F(Site01Nmer2, mir1_def) {
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
    EXPECT_EQ(5u, sites.get_length());

    test_sites2(sites, 0, "7mer", 0, 25, false);
    test_sites2(sites, 1, "7mer", 1, 25, false);
    test_sites2(sites, 2, "7mer", 2, 25, true);
    test_sites2(sites, 3, "7mer", 3, 25, true);
    test_sites2(sites, 4, "7mer", 4, 25, true);
}
}
