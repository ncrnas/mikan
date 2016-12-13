#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_rh2.hpp"

namespace {

    class Site03MM2 : public TestSiteRH2
    {
    protected:
        Site03MM2() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_03_mm_2.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_rh2/";

            mSeedDef1 = "6mGU1";
            mOverlapDef = "orig";
        }

        typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TIndexQGram TIdx;
        typedef rh2mfe::RH2Core<rh2mfe::TRNATYPE>::TFinder TFin;
        typedef rh2mfe::RH2SeedSites<rh2mfe::TRNATYPE> TSit;

    };

    TEST_F(Site03MM2, mir1_mm) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder, mrna_seqs);

        int ret_val = sites.find_seed_sites(mirna_seqs[1], mSeedDef1, mOverlapDef);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(4u, sites.get_length());

        test_sites2(sites, 0, "6mer", 0, 24, true);
        test_sites2(sites, 1, "6mer", 1, 24, true);
        test_sites2(sites, 2, "6mer", 2, 24, true);
        test_sites2(sites, 3, "6mer", 3, 24, true);
    }

}
