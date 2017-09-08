#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site05BT1 : public TestSiteMKENS {
protected:
    Site05BT1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_1.fasta";
    }

};

TEST_F(Site05BT1, mir124) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(0);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(27u, sites.get_length());

    test_sites(sites, 0, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 0, 30, true, 0);
    test_sites(sites, 1, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:BT,", 1, 30, true, 0);
    test_sites(sites, 2, "mr:8mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 2, 30, true, 0);
    test_sites(sites, 3, "mr:8mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:BT,", 3, 30, true, 0);
    test_sites(sites, 4, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 4, 30, true, 0);
    test_sites(sites, 5, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 5, 30, true, 0);
    test_sites(sites, 6, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 6, 30, true, 0);
    test_sites(sites, 7, "mr:8mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 7, 30, true, 0);
    test_sites(sites, 8, "mr:8mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:BT,", 8, 30, true, 0);
    test_sites(sites, 9, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:BT,", 9, 30, true, 0);
    test_sites(sites, 10, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT,", 10, 30, true, 0);
    test_sites(sites, 11, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 11, 30, true, 0);
    test_sites(sites, 12, "mr:8mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 12, 30, true, 0);
    test_sites(sites, 13, "mr:8mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:BT,", 13, 30, true, 0);
    test_sites(sites, 14, "mr:8mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:BT,", 14, 30, true, 0);
    test_sites(sites, 15, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT,", 15, 30, true, 0);
    test_sites(sites, 16, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 16, 30, true, 0);
    test_sites(sites, 17, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:NA,", 17, 30, true, 0);
    test_sites(sites, 18, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:BT,", 18, 30, true, 0);
    test_sites(sites, 19, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:BT,", 19, 30, true, 0);
    test_sites(sites, 20, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM,", 20, 29, true, 0);
    test_sites(sites, 21, "mr:7mer_BT,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA,", 20, 30, true, 0);
    test_sites(sites, 22, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM,", 21, 29, true, 0);
    test_sites(sites, 23, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:NA,", 21, 30, true, 0);
    test_sites(sites, 24, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:BT,", 22, 30, true, 0);
    test_sites(sites, 25, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:BT,", 23, 30, true, 0);
    test_sites(sites, 26, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:BT,", 24, 30, true, 0);

}

}
