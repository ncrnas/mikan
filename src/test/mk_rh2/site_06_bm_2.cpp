#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

class Site06BM2 : public TestSiteRH2 {
protected:
    Site06BM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_rh2/";

        mSeedDef1 = "6mer";
        mOverlapDef = "orig";
    }

    typedef rh2mfe::RH2Core<mikan::TRNATYPE>::TIndexQGram TIdx;
    typedef rh2mfe::RH2Core<mikan::TRNATYPE>::TFinder TFin;
    typedef rh2mfe::RH2SeedSites<mikan::TRNATYPE> TSit;

};

TEST_F(Site06BM2, mir1_bm_6mer) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef1, mOverlapDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(0u, sites.get_length());
}

TEST_F(Site06BM2, mir1_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef1 = "7mGU+";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef1, mOverlapDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(8u, sites.get_length());

    test_sites2(sites, 0, "7mer_GUT", 0, 15, true);
    test_sites2(sites, 1, "7mer_GUT", 1, 15, true);
    test_sites2(sites, 2, "7mer_GUT", 2, 15, true);
    test_sites2(sites, 3, "7mer_GUT", 3, 15, true);
    test_sites2(sites, 4, "7mer_GUT", 4, 15, true);
    test_sites2(sites, 5, "7mer_GUT", 5, 15, true);
    test_sites2(sites, 6, "7mer_GUT", 6, 15, true);
    test_sites2(sites, 7, "7mer_GUT", 7, 15, true);
}

}
