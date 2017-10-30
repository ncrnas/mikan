#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site06BM1 : public TestSiteMKENS {
protected:
    Site06BM1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_1.fasta";
    }

};

TEST_F(Site06BM1, mir124) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(0);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(8u, sites.get_length());

    test_sites(sites, 0, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 1, 31, true, 0);
    test_sites(sites, 1, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 3, 31, true, 0);
    test_sites(sites, 2, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 4, 32, true, 0);
    test_sites(sites, 3, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 5, 31, true, 0);
    test_sites(sites, 4, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 5, 32, true, 0);
    test_sites(sites, 5, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:LP", 6, 32, true, 0);
    test_sites(sites, 6, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 7, 31, true, 0);
    test_sites(sites, 7, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:NA", 7, 32, true, 0);

}

}
