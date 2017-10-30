#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site05BT2 : public TestSiteMKENS {
protected:
    Site05BT2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_05_bt_2.fasta";
    }

};

TEST_F(Site05BT2, mir1) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(1);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(55u, sites.get_length());

    test_sites(sites, 0, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 0, 21, true, 0);
    test_sites(sites, 1, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 0, 30, true, 0);
    test_sites(sites, 2, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 1, 21, true, 0);
    test_sites(sites, 3, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:BT", 1, 30, true, 0);
    test_sites(sites, 4, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 2, 21, true, 0);
    test_sites(sites, 5, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 2, 30, true, 0);
    test_sites(sites, 6, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 3, 21, true, 0);
    test_sites(sites, 7, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT", 3, 30, true, 0);
    test_sites(sites, 8, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 4, 21, true, 0);
    test_sites(sites, 9, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 4, 30, true, 0);
    test_sites(sites, 10, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 5, 21, true, 0);
    test_sites(sites, 11, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 5, 30, true, 0);
    test_sites(sites, 12, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 5, 31, true, 0);
    test_sites(sites, 13, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 6, 21, true, 0);
    test_sites(sites, 14, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 6, 30, true, 0);
    test_sites(sites, 15, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 7, 21, true, 0);
    test_sites(sites, 16, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 7, 30, true, 0);
    test_sites(sites, 17, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 8, 21, true, 0);
    test_sites(sites, 18, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT", 8, 30, true, 0);
    test_sites(sites, 19, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 9, 21, true, 0);
    test_sites(sites, 20, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT", 9, 30, true, 0);
    test_sites(sites, 21, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 10, 21, true, 0);
    test_sites(sites, 22, "mr:8mer_GUT,pt:8mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:BT", 10, 30, true, 0);
    test_sites(sites, 23, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 11, 21, true, 0);
    test_sites(sites, 24, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 11, 30, true, 0);
    test_sites(sites, 25, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 11, 32, true, 0);
    test_sites(sites, 26, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 12, 21, true, 0);
    test_sites(sites, 27, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 12, 30, true, 0);
    test_sites(sites, 28, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 12, 32, true, 0);
    test_sites(sites, 29, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 13, 21, true, 0);
    test_sites(sites, 30, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT", 13, 30, true, 0);
    test_sites(sites, 31, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 14, 21, true, 0);
    test_sites(sites, 32, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT", 14, 30, true, 0);
    test_sites(sites, 33, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 15, 21, true, 0);
    test_sites(sites, 34, "mr:8mer_BT,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:BT", 15, 30, true, 0);
    test_sites(sites, 35, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 16, 21, true, 0);
    test_sites(sites, 36, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 16, 30, true, 0);
    test_sites(sites, 37, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 17, 21, true, 0);
    test_sites(sites, 38, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:NA", 17, 30, true, 0);
    test_sites(sites, 39, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 18, 21, true, 0);
    test_sites(sites, 40, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:BT", 18, 30, true, 0);
    test_sites(sites, 41, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 19, 21, true, 0);
    test_sites(sites, 42, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:BT", 19, 30, true, 0);
    test_sites(sites, 43, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 20, 21, true, 0);
    test_sites(sites, 44, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 20, 29, true, 0);
    test_sites(sites, 45, "mr:7mer_BT,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 20, 30, true, 0);
    test_sites(sites, 46, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 21, 21, true, 0);
    test_sites(sites, 47, "mr:NA,pt:NA,rh:NA,tm:NA,ts:NA,sv:BM", 21, 29, true, 0);
    test_sites(sites, 48, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:NA", 21, 30, true, 0);
    test_sites(sites, 49, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 22, 21, true, 0);
    test_sites(sites, 50, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:BT", 22, 30, true, 0);
    test_sites(sites, 51, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 23, 21, true, 0);
    test_sites(sites, 52, "mr:8mer_BT,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:BT", 23, 30, true, 0);
    test_sites(sites, 53, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA", 24, 21, true, 0);
    test_sites(sites, 54, "mr:8mer_GUT,pt:8mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:BT", 24, 30, true, 0);

}

}
