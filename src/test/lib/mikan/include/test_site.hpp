#ifndef MIKAN_TEST_SITE_HPP_
#define MIKAN_TEST_SITE_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "test_main_io.hpp"

template<class TOpts, class TSeeds, class TSites>
class TestSite : public TestIOBase {
protected:

    void test_sites(const TSites &sites, int idx, const char *seq_type, unsigned mpos, unsigned spos,
                    bool effective, int mmpos) {
        seqan::String<unsigned> const &mrnapos = sites.get_mrna_pos();
        seqan::String<unsigned> const &sitepos = sites.get_site_pos();
        seqan::StringSet<seqan::CharString> const &seedtypes = sites.get_seed_types();
        seqan::String<int> const &mismatchpos = sites.get_mismatched_pos();

        EXPECT_EQ(mpos, mrnapos[idx]);
        EXPECT_EQ(spos, sitepos[idx]);
        EXPECT_STREQ(seq_type, (const char *) seqan::toCString(seedtypes[idx]));
        EXPECT_EQ(effective, sites.mEffectiveSites[idx]);
        EXPECT_EQ(mmpos, mismatchpos[idx]);

//        std::cout << "SS"  << idx <<  ", ";
//        std::cout << "CC" <<  (const char *)seqan::toCString(seedtypes[idx]) << "CC, ";
//        std::cout << mrnapos[idx] <<  ", ";
//        std::cout << sitepos[idx] <<  ", ";
//        std::cout << sites.mEffectiveSites[idx] <<  ", ";
//        std::cout <<  mismatchpos[idx] <<  ");";
//        std::cout << std::endl;
    }

    void test_sites2(const TSites &sites, int idx, const char *seq_type, unsigned mpos, unsigned spos,
                     bool effective) {
        seqan::String<unsigned> const &mrnapos = sites.get_mrna_pos();
        seqan::String<unsigned> const &sitepos = sites.get_site_pos();
        seqan::StringSet<seqan::CharString> const &seedtypes = sites.get_seed_types();

        EXPECT_EQ(mpos, mrnapos[idx]);
        EXPECT_EQ(spos, sitepos[idx]);
        EXPECT_STREQ(seq_type, (const char *) seqan::toCString(seedtypes[idx]));
        EXPECT_EQ(effective, sites.mEffectiveSites[idx]);
    }

    void test_sites3(const TSites &sites, int idx, unsigned mpos, unsigned spos) {
        seqan::String<unsigned> const &mrnapos = sites.get_mrna_pos();
        seqan::String<unsigned> const &sitepos = sites.get_site_pos();

        EXPECT_EQ(mpos, mrnapos[idx]);
        EXPECT_EQ(spos, sitepos[idx]);
    }

    seqan::StringSet<seqan::CharString> mSeedDef;
    seqan::CharString mSeedDef1;
    seqan::CharString mOverlapDef;

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef TSites TSit;
    typedef TSeeds TSeed;
    typedef TOpts TOp;

};

#endif //MIKAN_TEST_SITE_HPP_
