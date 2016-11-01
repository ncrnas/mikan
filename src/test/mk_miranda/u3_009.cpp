#include <iostream>
#include "gtest/gtest.h"
#include "test_io.hpp"
#include "mr3_core.hpp"

namespace {

    class U3009 : public TestIOMR3AS
    {
    protected:
        U3009() {
            IFNAME1 = (char *)"mir_001.fasta";
            IFNAME2 = (char *)"utr3_009.fasta";
            O1FNAME1 = (char *)"test_output1_site_9.txt";
            O1FNAME2 = (char *)"test_output1_mrna_9.txt";
            O2FNAME1 = (char *)"test_output2_site_9.txt";
            O2FNAME2 = (char *)"test_output2_mrna_9.txt";
            OMPATH = (char *)"mk_miranda/";
        }
    };

    TEST_F(U3009, comp_site) {
        run_main();

        gtest_compare_two_files(o1file1, o2file1);
    }

    TEST_F(U3009, comp_mrna) {
        run_main();

        gtest_compare_two_files(o1file2, o2file2);
    }
}