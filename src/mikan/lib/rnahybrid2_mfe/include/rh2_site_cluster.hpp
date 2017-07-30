#ifndef RH2_SITE_CLUSTER_HPP_
#define RH2_SITE_CLUSTER_HPP_

#include <set>                   // set
#include <map>                   // multimap
#include <utility>               // pair
#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_filter.hpp"    // MKSiteFilter
#include "mk_option.hpp"         // MKOptions
#include "rh2_score.hpp"         // RH2SiteScores
#include "rh2_seed_site.hpp"     // RH2SeedSites

namespace rh2mfe {

//
// RH2SiteCluster
//
class RH2SiteCluster {
public:
    // Define methods
    RH2SiteCluster() : mSiteCount(0) {}

    std::set<unsigned> &get_mrna_pos_set() { return mRNAPosSet; }

    std::multimap<unsigned, unsigned> &get_mrna_pos_map() { return mSiteMap; }

    int get_site_count() { return mSiteCount; }

    // Method prototype
    void clear_cluster();

    void cluster_site_pos(RH2SeedSites &pSeedSites, RH2SiteScores &pScores);

private:
    typedef std::pair<unsigned, unsigned> TPosPair;

    int mSiteCount;
    std::set<unsigned> mRNAPosSet;
    std::multimap<unsigned, unsigned> mSiteMap;

};

//
// Identify overlapped sites
//
class RH2SiteFilter : public mikan::MKSiteFilter {
public:
    // Define methods
    RH2SiteFilter(mikan::MKOptions const &opts) : MKSiteFilter(opts) {}

    // Method prototype
    void init_from_args();

private:
    float get_precedence(unsigned pSitePos, mikan::MKSeedSites &pSeedSites,
                         mikan::MKSiteScores &pSiteScores);

    void set_intervals(mikan::MKSeedSites &pSeedSites, mikan::MKSiteScores &pSiteScores,  unsigned pSiteIdx,
                       unsigned &pStartSearch, unsigned &pEndSearch, unsigned &pStartAdd, unsigned &pEndAdd,
                       bool &pSearchOverlap);

    seqan::CharString mOverlapMethod;

};

//
// Filter sites by the number of sites
//
class RH2TopNScore {
public:
    // Define methods
    RH2TopNScore() : mTopN(0) {}

    // Method prototype
    int filter_sites(RH2SeedSites &pSeedSites, RH2SiteScores &pScores, int pMaxHits);

    void clear_cluster();

private:
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<TItMap, TItMap> TItRetPair;
    typedef std::multimap<float, unsigned>::iterator TITStartScore;
    typedef std::pair<float, unsigned> TScorePair;

    RH2SiteCluster mSiteCluster;
    int mTopN;

private:
    void sort_sites_by_score(RH2SiteScores &pScores, int pPosIdx);

    void mark_non_topn_sites(RH2SiteScores &pScores, std::multimap<float, unsigned> &pSortedSites);


};

} // namespace rh2mfe

#endif /* RH2_SITE_CLUSTER_HPP_ */
