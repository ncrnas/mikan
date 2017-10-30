#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site01Nmer1 : public TestSiteMKENS {
protected:
    Site01Nmer1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_1.fasta";
    }
};

TEST_F(Site01Nmer1, mir124) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(0);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(53u, sites.get_length());

    test_sites(sites, 0, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 1, 7, true, 0);
    test_sites(sites, 1, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 2, 19, true, 0);
    test_sites(sites, 2, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 3, 20, true, 0);
    test_sites(sites, 3, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:6mer", 4, 21, true, 0);
    test_sites(sites, 4, "mr:7mer_MM,pt:6mer,rh:NA,tm:NA,ts:NA,sv:6mer", 5, 24, true, 0);
    test_sites(sites, 5, "mr:7mer_MM,pt:6mer,rh:NA,tm:NA,ts:NA,sv:6mer", 6, 25, true, 0);
    test_sites(sites, 6, "mr:7mer_MM,pt:6mer,rh:NA,tm:6mer,ts:NA,sv:6mer", 7, 26, true, 0);
    test_sites(sites, 7, "mr:7mer_MM,pt:6mer,rh:NA,tm:6mer,ts:NA,sv:6mer", 8, 38, true, 0);
    test_sites(sites, 8, "mr:7mer_MM,pt:6mer,rh:NA,tm:6mer,ts:NA,sv:6mer", 9, 39, true, 0);
    test_sites(sites, 9, "mr:7mer_MM,pt:NA,rh:NA,tm:6mer,ts:NA,sv:7mer-A1", 10, 40, true, 0);
    test_sites(sites, 10, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 12, 7, true, 0);
    test_sites(sites, 11, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 13, 19, true, 0);
    test_sites(sites, 12, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:7mer-A1,sv:NA", 14, 20, true, 0);
    test_sites(sites, 13, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:7mer-A1,sv:7mer-A1", 15, 21, true, 0);
    test_sites(sites, 14, "mr:7mer_MM,pt:6mer,rh:NA,tm:NA,ts:7mer-A1,sv:7mer-A1", 16, 24, true, 0);
    test_sites(sites, 15, "mr:7mer_MM,pt:6mer,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1", 17, 25, true, 0);
    test_sites(sites, 16, "mr:7mer_MM,pt:6mer,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1", 18, 26, true, 0);
    test_sites(sites, 17, "mr:7mer_MM,pt:6mer,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1", 19, 38, true, 0);
    test_sites(sites, 18, "mr:7mer_MM,pt:6mer,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1", 20, 39, true, 0);
    test_sites(sites, 19, "mr:7mer_MM,pt:NA,rh:NA,tm:6mer,ts:NA,sv:7mer-A1", 21, 40, true, 0);
    test_sites(sites, 20, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 22, 7, true, 0);
    test_sites(sites, 21, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 23, 19, true, 0);
    test_sites(sites, 22, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 24, 20, true, 0);
    test_sites(sites, 23, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:7mer-m8,sv:7mer-m8", 25, 21, true, 0);
    test_sites(sites, 24, "mr:7mer,pt:7mer,rh:7mer,tm:NA,ts:7mer-m8,sv:7mer-m8", 26, 24, true, 0);
    test_sites(sites, 25, "mr:7mer,pt:7mer,rh:7mer,tm:NA,ts:7mer-m8,sv:7mer-m8", 27, 25, true, 0);
    test_sites(sites, 26, "mr:7mer,pt:7mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 28, 26, true, 0);
    test_sites(sites, 27, "mr:7mer,pt:7mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 29, 27, true, 0);
    test_sites(sites, 28, "mr:7mer,pt:7mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 30, 38, true, 0);
    test_sites(sites, 29, "mr:7mer,pt:7mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 31, 39, true, 0);
    test_sites(sites, 30, "mr:7mer,pt:NA,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:8mer", 32, 40, true, 0);
    test_sites(sites, 31, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 33, 7, true, 0);
    test_sites(sites, 32, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 34, 19, true, 0);
    test_sites(sites, 33, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 35, 20, true, 0);
    test_sites(sites, 34, "mr:7mer,pt:NA,rh:7mer,tm:NA,ts:8mer,sv:8mer", 36, 21, true, 0);
    test_sites(sites, 35, "mr:7mer,pt:7mer,rh:7mer,tm:NA,ts:8mer,sv:8mer", 37, 24, true, 0);
    test_sites(sites, 36, "mr:7mer,pt:7mer,rh:7mer,tm:8mer,ts:8mer,sv:8mer", 38, 25, true, 0);
    test_sites(sites, 37, "mr:7mer,pt:7mer,rh:7mer,tm:8mer,ts:8mer,sv:8mer", 39, 26, true, 0);
    test_sites(sites, 38, "mr:7mer,pt:7mer,rh:7mer,tm:8mer,ts:8mer,sv:8mer", 40, 27, true, 0);
    test_sites(sites, 39, "mr:7mer,pt:7mer,rh:7mer,tm:8mer,ts:8mer,sv:8mer", 41, 38, true, 0);
    test_sites(sites, 40, "mr:7mer,pt:7mer,rh:7mer,tm:8mer,ts:8mer,sv:8mer", 42, 39, true, 0);
    test_sites(sites, 41, "mr:7mer,pt:NA,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:8mer", 43, 40, true, 0);
    test_sites(sites, 42, "mr:8mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 44, 8, true, 0);
    test_sites(sites, 43, "mr:8mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 45, 19, true, 0);
    test_sites(sites, 44, "mr:8mer,pt:NA,rh:7mer,tm:NA,ts:NA,sv:NA", 46, 20, true, 0);
    test_sites(sites, 45, "mr:8mer,pt:NA,rh:7mer,tm:NA,ts:7mer-m8,sv:7mer-m8", 47, 21, true, 0);
    test_sites(sites, 46, "mr:8mer,pt:8mer,rh:7mer,tm:NA,ts:7mer-m8,sv:7mer-m8", 48, 24, true, 0);
    test_sites(sites, 47, "mr:8mer,pt:8mer,rh:7mer,tm:NA,ts:7mer-m8,sv:7mer-m8", 49, 25, true, 0);
    test_sites(sites, 48, "mr:8mer,pt:8mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 50, 26, true, 0);
    test_sites(sites, 49, "mr:8mer,pt:8mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 51, 27, true, 0);
    test_sites(sites, 50, "mr:8mer,pt:8mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 52, 38, true, 0);
    test_sites(sites, 51, "mr:8mer,pt:8mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8", 53, 39, true, 0);
    test_sites(sites, 52, "mr:8mer,pt:NA,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:8mer", 54, 40, true, 0);

}

}
