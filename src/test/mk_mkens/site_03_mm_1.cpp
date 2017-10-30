#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site03MM1 : public TestSiteMKENS {
protected:
    Site03MM1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_03_mm_1.fasta";
    }

};

TEST_F(Site03MM1, mir124) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(0);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(28u, sites.get_length());

    test_sites(sites, 0, "mr:7mer_MM,pt:6mer,rh:NA,tm:6mer,ts:NA,sv:6mer", 0, 30, true, 0);
    test_sites(sites, 1, "mr:7mer_MM,pt:6mer,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1", 1, 30, true, 0);
    test_sites(sites, 2, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:6mer", 2, 30, true, 0);
    test_sites(sites, 3, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1", 3, 30, true, 0);
    test_sites(sites, 4, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 6, 30, true, 0);
    test_sites(sites, 5, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 7, 30, true, 0);
    test_sites(sites, 6, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:NA", 8, 30, true, 0);
    test_sites(sites, 7, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:LP", 9, 30, true, 0);
    test_sites(sites, 8, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 12, 30, true, 0);
    test_sites(sites, 9, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 13, 30, true, 0);
    test_sites(sites, 10, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:NA", 14, 30, true, 0);
    test_sites(sites, 11, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:LP", 15, 30, true, 0);
    test_sites(sites, 12, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 18, 30, true, 0);
    test_sites(sites, 13, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 19, 30, true, 0);
    test_sites(sites, 14, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:NA", 20, 30, true, 0);
    test_sites(sites, 15, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:LP", 21, 30, true, 0);
    test_sites(sites, 16, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 24, 30, true, 0);
    test_sites(sites, 17, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 25, 30, true, 0);
    test_sites(sites, 18, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:NA", 26, 30, true, 0);
    test_sites(sites, 19, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:LP", 27, 30, true, 0);
    test_sites(sites, 20, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA", 30, 30, true, 0);
    test_sites(sites, 21, "mr:7mer_MM,pt:NA,rh:NA,tm:NA,ts:NA,sv:LP", 31, 30, true, 0);
    test_sites(sites, 22, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:NA", 32, 30, true, 0);
    test_sites(sites, 23, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:NA,ts:NA,sv:LP", 33, 30, true, 0);
    test_sites(sites, 24, "mr:7mer_MM,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA", 36, 30, true, 0);
    test_sites(sites, 25, "mr:7mer_MM,pt:NA,rh:NA,tm:6mer,ts:NA,sv:LP", 37, 30, true, 0);
    test_sites(sites, 26, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:NA", 38, 30, true, 0);
    test_sites(sites, 27, "mr:8mer_MM,pt:8mer_MM,rh:NA,tm:6mer,ts:NA,sv:LP", 39, 30, true, 0);

}
}
