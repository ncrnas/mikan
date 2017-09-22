#include <iostream>
#include "gtest/gtest.h"
#include "test_tm1.hpp"

namespace {

class OM001U001 : public TestIOBase {
protected:
    OM001U001() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_tm1/";
    }
};

TEST_F(OM001U001, comp_site) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    //TODO: Need to recheck hsa-miR-1 NM_000165 687 693 6mer - only mikan finds this site
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM001U001, comp_mrna) {
    (void) mikan::MKMain<tm1p::TM1Options, tm1p::TM1Core>(argc, (const char **) argv);
    //TODO: Need to recheck hsa-miR-1 NM_000165 687 693 6mer - only mikan finds this site
    gtest_compare_two_files(o1file2, o2file2);
}
}