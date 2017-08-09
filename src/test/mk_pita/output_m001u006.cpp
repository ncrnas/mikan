#include <iostream>
#include "gtest/gtest.h"
#include "test_pita.hpp"

namespace {

class OM001U006 : public TestIOPITA {
protected:
    OM001U006() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_006.fasta";
        O1FNAME1 = (char *) "test_output1_site_6.txt";
        O1FNAME2 = (char *) "test_output1_mrna_6.txt";
        O2FNAME1 = (char *) "test_output2_site_6.txt";
        O2FNAME2 = (char *) "test_output2_mrna_6.txt";
        OMPATH = (char *) "mk_pita/";
    }
};

TEST_F(OM001U006, comp_site) {
    (void) mikan::MKCoreMain<ptddg::PITAOptions, ptddg::PITACore >(argc, (const char **) argv);
    gtest_compare_two_files(o1file1, o2file1);
}

TEST_F(OM001U006, comp_mrna) {
    (void) mikan::MKCoreMain<ptddg::PITAOptions, ptddg::PITACore >(argc, (const char **) argv);
    gtest_compare_two_files(o1file2, o2file2);
}
}