#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

    class Site02GU1 : public TestSiteTS5
    {
    protected:
        Site02GU1() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"ts_02_gu_1.fasta";
            O1FNAME1 = (char *)"test_output1_site_1.txt";
            O1FNAME2 = (char *)"test_output1_mrna_1.txt";
            O2FNAME1 = (char *)"test_output2_site_1.txt";
            O2FNAME2 = (char *)"test_output2_mrna_1.txt";
            OMPATH = (char *)"mk_ts5/";
        }

        typedef ts5cs::TS5Core<ts5cs::TRNATYPE>::TIndexQGram TIdx;
        typedef ts5cs::TS5Core<ts5cs::TRNATYPE>::TFinder TFin;
        typedef ts5cs::TS5SeedSites<ts5cs::TRNATYPE> TSit;

    };

    TEST_F(Site02GU1, mir124_gu) {
        read_files(false);
        set_seqs();
        TIdx index(mrna_seqs);
        TFin finder(index);
        TSit sites(index, finder);

        int ret_val = sites.find_seed_sites(mirna_seqs[0]);
        EXPECT_EQ(0, ret_val);
        EXPECT_EQ(1u, sites.get_length());

        test_sites3(sites, 0, 0, 24);
    }
}