#ifndef TM1_SITE_CLUSTER_HPP_
#define TM1_SITE_CLUSTER_HPP_

#include <tm1_seed_site.hpp>      // TM1SeedSites
#include <set>                    // set
#include <map>                    // multimap
#include <utility>                // pair

namespace tm1p {

//
// Provide basic site cluster
//
template<class TRNAString>
class TM1SiteCluster {
public:
    // Define methods
    TM1SiteCluster() : mSiteCount(0) {}

    std::set<unsigned> &get_mrna_pos_set() { return mRNAPosSet; }

    std::multimap<unsigned, unsigned> &get_mrna_pos_map() { return mSiteMap; }

    int get_site_count() { return mSiteCount; }

    // Method prototype
    void clear_cluster();

    void cluster_site_pos(TM1SeedSites <TRNAString> &pSeedSites);

private:
    typedef std::pair<unsigned, unsigned> TPosPair;

    int mSiteCount;
    std::set<unsigned> mRNAPosSet;
    std::multimap<unsigned, unsigned> mSiteMap;

};

//
// Sort sites by position
//
template<class TRNAString>
class TM1SortedSitePos {
public:
    // Define methods
    TM1SortedSitePos() {}

    // Method prototype
    int generate_sorted_mrna_pos(TM1SeedSites <TRNAString> &pSeedSites, bool pRemoveOvelaps);

    const seqan::StringSet<seqan::String<unsigned> > &get_sorted_mrna_pos() { return mSortedSites; }

    const seqan::String<unsigned> &get_mrna_ids() { return mMRNAIDs; }

    void clear_site_pos();

private:
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMMap;
    typedef std::map<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<TItMMap, TItMMap> TItRetPair;
    typedef std::multimap<unsigned, unsigned>::iterator TITStartPos;
    typedef std::pair<unsigned, unsigned> TPosPair;

    seqan::StringSet<seqan::String<unsigned> > mSortedSites;
    seqan::String<unsigned> mMRNAIDs;
    TM1SiteCluster<TRNAString> mSiteCluster;

private:
    void remove_overlapped_sites(TM1SeedSites <TRNAString> &pSeedSites);

    void sort_by_seed_types(TM1SeedSites <TRNAString> &pSeedSites, TItRetPair &pGroupedSites,
                            std::map<unsigned, unsigned> &pSortedSites);

    void sort_by_pos7(TM1SeedSites <TRNAString> &pSeedSites, TItRetPair &pGroupedSites,
                      std::map<unsigned, unsigned> &pSortedSites);

};

} // namespace tm1p

#endif /* TM1_SITE_CLUSTER_HPP_ */
