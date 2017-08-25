#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class SeedAll : public TestSeedTS5 {
protected:
    SeedAll() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "utr3_001.fasta";
    }
};

TEST_F(SeedAll, mir124_def) {
    create_seed_seqs(0);

    test_seed("AAGGCA", 0, "6mer", true, 0);

}

TEST_F(SeedAll, mir1_def) {
    create_seed_seqs(1);

    test_seed("GGAAUG", 0, "6mer", true, 0);

}

}
