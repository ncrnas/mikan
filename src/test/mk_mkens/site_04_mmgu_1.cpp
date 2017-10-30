#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site04MMGU1 : public TestSiteMKENS {
protected:
    Site04MMGU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_04_mmgu_1.fasta";
    }

};

TEST_F(Site04MMGU1, mir124) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(0);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(66u, sites.get_length());

    test_sites(sites, 0, "mr:7mer_MM,pt:6mer,rh:NA,tm:6mer,ts:NA,sv:6mer", 0, 30, true, 0);
    test_sites(sites, 1, "mr:7mer_MM,pt:6mer,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1", 1, 30, true, 0);
    test_sites(sites, 2, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 2, 30, true, 0);
    test_sites(sites, 3, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:7mer-A1,ts:NA,sv:NA", 3, 30, true, 0);
    test_sites(sites, 4, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:NA", 4, 30, true, 0);
    test_sites(sites, 5, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:8mer,ts:NA,sv:GUM", 5, 30, true, 0);
    test_sites(sites, 6, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 6, 30, true, 0);
    test_sites(sites, 7, "mr:7mer_MMGU,pt:NA,rh:NA,tm:7mer-A1,ts:NA,sv:NA", 7, 30, true, 0);
    test_sites(sites, 8, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 8, 30, true, 0);
    test_sites(sites, 9, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:7mer-A1,ts:NA,sv:NA", 9, 30, true, 0);
    test_sites(sites, 10, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:NA", 10, 30, true, 0);
    test_sites(sites, 11, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:8mer,ts:NA,sv:GUM", 11, 30, true, 0);
    test_sites(sites, 12, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 12, 30, true, 0);
    test_sites(sites, 13, "mr:7mer_MMGU,pt:NA,rh:NA,tm:7mer-A1,ts:NA,sv:NA", 13, 30, true, 0);
    test_sites(sites, 14, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 14, 30, true, 0);
    test_sites(sites, 15, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 15, 30, true, 0);
    test_sites(sites, 16, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 16, 30, true, 0);
    test_sites(sites, 17, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 17, 30, true, 0);
    test_sites(sites, 18, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 18, 30, true, 0);
    test_sites(sites, 19, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 19, 30, true, 0);
    test_sites(sites, 20, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 20, 30, true, 0);
    test_sites(sites, 21, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 21, 30, true, 0);
    test_sites(sites, 22, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 22, 30, true, 0);
    test_sites(sites, 23, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 23, 30, true, 0);
    test_sites(sites, 24, "mr:7mer_MM,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 24, 30, true, 0);
    test_sites(sites, 25, "mr:7mer_MM,pt:NA,rh:NA,tm:6mer,ts:NA,sv:LP", 25, 30, true, 0);
    test_sites(sites, 26, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 26, 30, true, 0);
    test_sites(sites, 27, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 27, 30, true, 0);
    test_sites(sites, 28, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 28, 30, true, 0);
    test_sites(sites, 29, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 29, 30, true, 0);
    test_sites(sites, 30, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 30, 30, true, 0);
    test_sites(sites, 31, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 31, 30, true, 0);
    test_sites(sites, 32, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 32, 30, true, 0);
    test_sites(sites, 33, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 33, 30, true, 0);
    test_sites(sites, 34, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 34, 30, true, 0);
    test_sites(sites, 35, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 35, 30, true, 0);
    test_sites(sites, 36, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 36, 30, true, 0);
    test_sites(sites, 37, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 37, 30, true, 0);
    test_sites(sites, 38, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 38, 30, true, 0);
    test_sites(sites, 39, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 39, 30, true, 0);
    test_sites(sites, 40, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 40, 30, true, 0);
    test_sites(sites, 41, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 41, 30, true, 0);
    test_sites(sites, 42, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 42, 30, true, 0);
    test_sites(sites, 43, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 43, 30, true, 0);
    test_sites(sites, 44, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 44, 30, true, 0);
    test_sites(sites, 45, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 45, 30, true, 0);
    test_sites(sites, 46, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 46, 30, true, 0);
    test_sites(sites, 47, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 47, 30, true, 0);
    test_sites(sites, 48, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 48, 30, true, 0);
    test_sites(sites, 49, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 49, 30, true, 0);
    test_sites(sites, 50, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 50, 30, true, 0);
    test_sites(sites, 51, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 51, 30, true, 0);
    test_sites(sites, 52, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 52, 30, true, 0);
    test_sites(sites, 53, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 53, 30, true, 0);
    test_sites(sites, 54, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 54, 30, true, 0);
    test_sites(sites, 55, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 55, 30, true, 0);
    test_sites(sites, 56, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 56, 30, true, 0);
    test_sites(sites, 57, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 57, 30, true, 0);
    test_sites(sites, 58, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 58, 30, true, 0);
    test_sites(sites, 59, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:NA,ts:NA,sv:NA", 59, 30, true, 0);
    test_sites(sites, 60, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 60, 30, true, 0);
    test_sites(sites, 61, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 61, 30, true, 0);
    test_sites(sites, 62, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 62, 30, true, 0);
    test_sites(sites, 63, "mr:8mer_MMGU,pt:8mer_MMGU,rh:NA,tm:6mer,ts:NA,sv:NA", 63, 30, true, 0);
    test_sites(sites, 64, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 64, 30, true, 0);
    test_sites(sites, 65, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 65, 30, true, 0);

}
}
