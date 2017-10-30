#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site06BM2 : public TestSiteMKENS {
protected:
    Site06BM2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_06_bm_2.fasta";
    }

};

TEST_F(Site06BM2, mir1) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(1);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(14u, sites.get_length());

    test_sites(sites, 0, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 0, 21, true, 0);
    test_sites(sites, 1, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 1, 21, true, 0);
    test_sites(sites, 2, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 1, 31, true, 0);
    test_sites(sites, 3, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 2, 21, true, 0);
    test_sites(sites, 4, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 3, 21, true, 0);
    test_sites(sites, 5, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 3, 31, true, 0);
    test_sites(sites, 6, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 4, 21, true, 0);
    test_sites(sites, 7, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 5, 21, true, 0);
    test_sites(sites, 8, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 5, 31, true, 0);
    test_sites(sites, 9, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 6, 21, true, 0);
    test_sites(sites, 10, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:LP", 6, 32, true, 0);
    test_sites(sites, 11, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 7, 21, true, 0);
    test_sites(sites, 12, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 7, 31, true, 0);
    test_sites(sites, 13, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:NA", 7, 32, true, 0);
}
}
