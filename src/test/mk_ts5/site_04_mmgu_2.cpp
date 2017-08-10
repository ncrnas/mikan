#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_ts5.hpp"

namespace {

class Site04MMGU2 : public TestSiteTS5 {
protected:
    Site04MMGU2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_2.fasta";
        O1FNAME1 = (char *) "test_output1_site_1.txt";
        O1FNAME2 = (char *) "test_output1_mrna_1.txt";
        O2FNAME1 = (char *) "test_output2_site_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        O2FNAME2 = (char *) "test_output2_mrna_1.txt";
        OMPATH = (char *) "mk_ts5/";
    }

};

TEST_F(Site04MMGU2, mir1_mmgu) {
    read_files();
    set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    TSit sites(index, finder, mrna_seqs);

    TOp ops;
    TSeed seedSeqs(ops);
    seedSeqs.set_seed_type_def(mSeedDef);
    seedSeqs.set_flags();;
    seedSeqs.create_seed_seqs(mirna_seqs[1]);

    int ret_val = sites.find_seed_sites(seedSeqs);
    EXPECT_EQ(0, ret_val);
    EXPECT_EQ(1u, sites.get_length());

    test_sites2(sites, 0, "7mer-A1", 1, 24, true);
}

}
