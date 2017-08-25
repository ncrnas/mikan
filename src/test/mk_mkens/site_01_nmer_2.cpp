#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site01Nmer2 : public TestSiteMKENS {
protected:
    Site01Nmer2() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_01_nmer_2.fasta";
    }
};

TEST_F(Site01Nmer2, mir1) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(1);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(14u, sites.get_length());

    test_sites(sites, 0, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 0, 19, true, 0);
    test_sites(sites, 1, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 0, 30, true, 0);
    test_sites(sites, 2, "mr:7mer_MM,pt:6mer,rh:NA,tm:6mer,ts:NA,sv:6mer,", 0, 31, true, 0);
    test_sites(sites, 3, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 1, 19, true, 0);
    test_sites(sites, 4, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 1, 30, true, 0);
    test_sites(sites, 5, "mr:7mer_MM,pt:6mer,rh:NA,tm:7mer-A1,ts:7mer-A1,sv:7mer-A1,", 1, 31, true, 0);
    test_sites(sites, 6, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 2, 19, true, 0);
    test_sites(sites, 7, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 2, 30, true, 0);
    test_sites(sites, 8, "mr:7mer,pt:7mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8,", 2, 31, true, 0);
    test_sites(sites, 9, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 3, 19, true, 0);
    test_sites(sites, 10, "mr:7mer_BT,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 3, 30, true, 0);
    test_sites(sites, 11, "mr:7mer,pt:7mer,rh:7mer,tm:8mer,ts:8mer,sv:8mer,", 3, 31, true, 0);
    test_sites(sites, 12, "mr:7mer_MMGU,pt:NA,rh:NA,tm:NA,ts:NA,sv:NA,", 4, 19, true, 0);
    test_sites(sites, 13, "mr:8mer,pt:8mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8,", 4, 31, true, 0);

}

}
