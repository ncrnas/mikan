#ifndef MIKAN_TEST_SITE_HPP_
#define MIKAN_TEST_SITE_HPP_

#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "test_main_io.hpp"

template<class TOptions, class TSeedSeqs, class TSeedSites>
class TestSite : public TestIOBase {
public:

    TestSite() : mSeedSeqs(mOpts) {}

protected:

    void test_sites(const TSeedSites &sites, unsigned idx, const char *seq_type, unsigned mpos, unsigned spos,
                    bool effective, int mmpos) {
        seqan::String<unsigned> const &mrnapos = sites.get_mrna_pos();
        seqan::String<unsigned> const &sitepos = sites.get_site_pos();
        seqan::StringSet<seqan::CharString> const &seedtypes = sites.get_seed_types();
        seqan::String<int> const &mismatchpos = sites.get_mismatched_pos();

        if (idx < length(mrnapos)) {
            EXPECT_EQ(mpos, mrnapos[idx]);
        }

        if (idx < length(sitepos)) {
            EXPECT_EQ(spos, sitepos[idx]);
        }

        if (idx < length(seedtypes)) {
            EXPECT_STREQ(seq_type, (const char *) seqan::toCString(seedtypes[idx]));
        }

        if (idx < length(sites.mEffectiveSites)) {
            EXPECT_EQ(effective, sites.mEffectiveSites[idx]);
        }

        if (idx < length(mismatchpos)) {
            EXPECT_EQ(mmpos, mismatchpos[idx]);
        }

    }

    void create_seed_seqs(unsigned pIdx) {
        read_and_set_seqs();

        mSeedSeqs.set_seed_type_def(mSeedDef);
        mSeedSeqs.set_flags();
        int ret = mSeedSeqs.create_seed_seqs(mirna_seqs[pIdx]);
        EXPECT_EQ(0, ret);
    }

    void find_seed_sites(TSeedSites &sites) {
        int ret = sites.find_seed_sites(mSeedSeqs);
        EXPECT_EQ(0, ret);
    }

    typedef mikan::TIndexQGram TIdx;
    typedef mikan::TFinder TFin;
    typedef TSeedSites TSit;

    TOptions mOpts;
    TSeedSeqs mSeedSeqs;
    seqan::StringSet<seqan::CharString> mSeedDef;
    seqan::CharString mOverlapDef;

};

#endif //MIKAN_TEST_SITE_HPP_
