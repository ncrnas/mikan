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

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef ts5cs::TS5SeedSites TSit;
    typedef ts5cs::TS5SeedSeqs TSeed;
    
};

TEST_F(Site01Nmer1, mir124) {
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
    EXPECT_EQ(31u, sites.get_length());

    test_sites2(sites, 1, "7mer-A1", 15, 15, true);
    test_sites2(sites, 2, "7mer-A1", 16, 18, true);
    test_sites2(sites, 3, "7mer-A1", 17, 19, true);
    test_sites2(sites, 4, "7mer-A1", 18, 20, true);
    test_sites2(sites, 5, "7mer-A1", 19, 32, true);
    test_sites2(sites, 6, "7mer-A1", 20, 33, true);
    test_sites2(sites, 7, "7mer-m8", 25, 15, true);
    test_sites2(sites, 8, "7mer-m8", 26, 18, true);
    test_sites2(sites, 9, "7mer-m8", 27, 19, true);
    test_sites2(sites, 10, "7mer-m8", 28, 20, true);
    test_sites2(sites, 11, "7mer-m8", 29, 21, true);
    test_sites2(sites, 12, "7mer-m8", 30, 32, true);
    test_sites2(sites, 13, "7mer-m8", 31, 33, true);
    test_sites2(sites, 14, "7mer-m8", 32, 34, true);
    test_sites2(sites, 15, "8mer", 36, 15, true);
    test_sites2(sites, 16, "8mer", 37, 18, true);
    test_sites2(sites, 17, "8mer", 38, 19, true);
    test_sites2(sites, 18, "8mer", 39, 20, true);
    test_sites2(sites, 19, "8mer", 40, 21, true);
    test_sites2(sites, 20, "8mer", 41, 32, true);
    test_sites2(sites, 21, "8mer", 42, 33, true);
    test_sites2(sites, 22, "7mer-m8", 43, 34, true);
    test_sites2(sites, 23, "7mer-m8", 47, 15, true);
    test_sites2(sites, 24, "7mer-m8", 48, 18, true);
    test_sites2(sites, 25, "7mer-m8", 49, 19, true);
    test_sites2(sites, 26, "7mer-m8", 50, 20, true);
    test_sites2(sites, 27, "7mer-m8", 51, 21, true);
    test_sites2(sites, 28, "7mer-m8", 52, 32, true);
    test_sites2(sites, 29, "7mer-m8", 53, 33, true);
    test_sites2(sites, 30, "7mer-m8", 54, 34, true);
}

}
