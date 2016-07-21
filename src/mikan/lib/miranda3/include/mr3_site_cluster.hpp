#ifndef MR3_SITE_CLUSTER_HPP_
#define MR3_SITE_CLUSTER_HPP_

#include <mikan/lib/miranda3/include/mr3_score.hpp>          // MR3DDGScores
#include <mikan/lib/miranda3/include/mr3_seed_site.hpp>      // MR3SeedSites
#include <set>                    // set
#include <map>                    // multimap
#include <utility>                // pair

namespace mr3as{

template <class TRNAString>

//
// Cluster target sites
//
class MR3SiteCluster
{
public:
    // Define methods
    MR3SiteCluster() : mSiteCount(0) {}
    std::set<unsigned>& get_mrna_pos_set() {return mRNAPosSet;}
    std::multimap<unsigned, unsigned>& get_mrna_pos_map() {return mSiteMap;}
    int get_site_count() {return mSiteCount;}

    // Method prototype
    void clear_cluster();
    void cluster_site_pos(MR3SeedSites<TRNAString> &pSeedSites);

private:
    typedef std::pair<unsigned, unsigned> TPosPair;

    int mSiteCount;
    std::set<unsigned> mRNAPosSet;
    std::multimap<unsigned, unsigned> mSiteMap;

};

//
// Identify overlapped sites
//
template <class TRNAString>
class MR3Overlap
{
public:
    // Define methods
    MR3Overlap() {}

    // Method prototype
    int filter_overlapped_sites(MR3SeedSites<TRNAString> &pSeedSites, int pGapLen);
    int filter_overlapped_sites_by_scores(MR3SeedSites<TRNAString> &pSeedSites,
            MR3SiteScores<TRNAString> &pSiteScores, int pGapLen);
    void clear_cluster();

private:
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<TItMap, TItMap> TItRetPair;
    typedef std::multimap<unsigned, unsigned>::iterator TITStartPos;
    typedef std::pair<unsigned, unsigned> TPosPair;

    MR3SiteCluster<TRNAString> mSiteCluster;

private:
    int make_overlapped_pairs(MR3SeedSites<TRNAString> &pSeedSites, int pGapLen,
            seqan::StringSet<seqan::String<unsigned> >& pPairs);
    void mark_overlapped_sites(MR3SeedSites<TRNAString> &pSeedSites, int pPrevIdx, int pCurIdx);
    unsigned get_seedtype_precedence(const seqan::CharString &pSeedType);
};

//
// Sort sites by position
//
template <class TRNAString>
class MR3SortedSitePos
{
public:
    // Define methods
    MR3SortedSitePos() {}

    // Method prototype
    int generate_sorted_mrna_pos(MR3SeedSites<TRNAString> &pSeedSites);
    const seqan::String<unsigned>& get_sorted_mrna_pos() {return mSortedSitePos;}
    void clear_site_pos();

private:
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<TItMap, TItMap> TItRetPair;
    typedef std::multimap<unsigned, unsigned>::iterator TITStartPos;
    typedef std::pair<unsigned, unsigned> TPosPair;

    seqan::String<unsigned> mSortedSitePos;
    MR3SiteCluster<TRNAString> mSiteCluster;

};

} // namespace mr3as

#endif /* MR3_SITE_CLUSTER_HPP_ */
