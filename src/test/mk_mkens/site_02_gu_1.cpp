#include<string>
#include <seqan/sequence.h>
#include "gtest/gtest.h"
#include "test_mkens.hpp"

namespace {

class Site02GU1 : public TestSiteMKENS {
protected:
    Site02GU1() {
        IFNAME1 = (char *) "mir_001.fasta";
        IFNAME2 = (char *) "ts_02_gu_1.fasta";
    }
};

TEST_F(Site02GU1, mir124) {
    read_and_set_seqs();
    TIdx index(mrna_seqs);
    TFin finder(index);
    mkens::MKECore mkCore(mOpts, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs, index, finder);
    mkCore.find_seed_sites(0);
    TSit &sites = mkCore.get_seed_sites();

    EXPECT_EQ(15u, sites.get_length());

    test_sites(sites, 0, "mr:7mer,pt:7mer,rh:7mer,tm:7mer-m8,ts:7mer-m8,sv:7mer-m8,", 0, 30, true, 0);
    test_sites(sites, 1, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA,", 1, 30, true, 0);
    test_sites(sites, 2, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:NA,", 2, 30, true, 0);
    test_sites(sites, 3, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:8mer,ts:NA,sv:GUM,", 3, 30, true, 0);
    test_sites(sites, 4, "mr:8mer_GUT,pt:8mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:NA,", 4, 30, true, 0);
    test_sites(sites, 5, "mr:8mer_GUT,pt:8mer_GUT,rh:7mer_GUT,tm:8mer,ts:NA,sv:GUM,", 5, 30, true, 0);
    test_sites(sites, 6, "mr:7mer_MMGU,pt:NA,rh:NA,tm:6mer,ts:NA,sv:NA,", 6, 30, true, 0);
    test_sites(sites, 7, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:NA,", 7, 30, true, 0);
    test_sites(sites, 8, "mr:7mer_GUT,pt:7mer_GUT,rh:7mer_GUT,tm:8mer,ts:NA,sv:GUM,", 8, 30, true, 0);
    test_sites(sites, 9, "mr:8mer_GUT,pt:8mer_GUT,rh:7mer_GUT,tm:7mer-m8,ts:NA,sv:NA,", 9, 30, true, 0);
    test_sites(sites, 10, "mr:8mer_GUT,pt:8mer_GUT,rh:7mer_GUT,tm:8mer,ts:NA,sv:GUM,", 10, 30, true, 0);
    test_sites(sites, 11, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA,", 12, 30, true, 0);
    test_sites(sites, 12, "mr:7mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA,", 13, 30, true, 0);
    test_sites(sites, 13, "mr:8mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA,", 14, 30, true, 0);
    test_sites(sites, 14, "mr:8mer_GU+,pt:NA,rh:7mer_GUT,tm:NA,ts:NA,sv:NA,", 15, 30, true, 0);

}

}
