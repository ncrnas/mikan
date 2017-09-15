#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class Site04MMGU1 : public TestSiteTS5 {
protected:
    Site04MMGU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_1.fasta";
    }

};

TEST_F(Site04MMGU1, mir124_mmgu) {
    create_seed_seqs(0);
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);
    find_seed_sites(sites);

    EXPECT_EQ(1u, sites.get_length());

    test_sites(sites, 0, "7mer-A1", 1, 24, true, 0);
}
}
