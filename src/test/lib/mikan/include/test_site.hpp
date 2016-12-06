#include<string>
#include "gtest/gtest.h"
#include "get_data_path.hpp"
#include "mikan_utils.hpp"
#include "mr3_core.hpp"

template <class TSites, class TTestIO>
class TestSite : public TTestIO
{
protected:

    void test_sites(const TSites &sites, int idx, const char *seq_type, unsigned mpos, unsigned spos,
                    bool effective, int mmpos) {
        seqan::String<unsigned> const& mrnapos = sites.get_mrna_pos();
        seqan::String<unsigned> const& sitepos = sites.get_site_pos();
        seqan::StringSet<seqan::CharString> const& seedtypes =  sites.get_seed_types();
        seqan::String<int> const& mismatchpos =  sites.get_mismatched_pos();

        EXPECT_EQ(mpos, mrnapos[idx]);
        EXPECT_EQ(spos, sitepos[idx]);
        EXPECT_STREQ(seq_type, (const char *)seqan::toCString(seedtypes[idx]));
        EXPECT_EQ(effective, sites.mEffectiveSites[idx]);
        EXPECT_EQ(mmpos, mismatchpos[idx]);
    }

    seqan::StringSet<seqan::CharString> mSeedDef;
};

typedef TestSite<mr3as::MR3SeedSites<mr3as::TRNATYPE>, TestIOMR3AS> TestSiteMR3AS;
