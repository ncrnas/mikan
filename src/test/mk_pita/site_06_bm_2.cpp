#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class Site06BM2 : public TestSitePITA {
protected:
    Site06BM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_2.fasta";
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
        mSeedDef[4] = "1:1";
        mSeedDef[5] = "1";
    }

    typedef ptddg::PITACore<ptddg::TRNATYPE>::TIndexQGram TIdx;
    typedef ptddg::PITACore<ptddg::TRNATYPE>::TFinder TFin;
    typedef ptddg::PITASeedSites<ptddg::TRNATYPE> TSit;

};

TEST_F(Site06BM2, mir1_bm) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(4u, sites.get_length());

    test_sites(sites, 0, "MM", 0, 25, false, 0);
    test_sites(sites, 1, "MM", 1, 25, false, 0);
    test_sites(sites, 2, "8mer_MM", 7, 26, true, 5);
    test_sites(sites, 3, "8mer_MM", 6, 26, true, 5);
}

TEST_F(Site06BM2, mir1_def) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    mSeedDef[3] = "1";
    mSeedDef[4] = "0:1";
    mSeedDef[5] = "1";
    int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(4u, sites.get_length());

    test_sites(sites, 0, "MM", 0, 25, false, 0);
    test_sites(sites, 1, "MM", 1, 25, false, 0);
    test_sites(sites, 2, "8mer_MM", 7, 26, true, 5);
    test_sites(sites, 3, "8mer_MM", 6, 26, true, 5);
}

}
