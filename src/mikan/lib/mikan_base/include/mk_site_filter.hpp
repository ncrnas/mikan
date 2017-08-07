#ifndef MK_SITE_FILTER_HPP_
#define MK_SITE_FILTER_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_option.hpp"            // MKOptions
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_site_score.hpp"        // MKSiteScores
#include "mk_rna_sites.hpp"         // MKRMAWithSites

namespace mikan {

//
// Filter overlapped sites
//
class MKSiteFilter {
public:
    // Define methods
    explicit MKSiteFilter(mikan::MKOptions const &opts) : mOpts(opts) {
        set_overlap_len(0);
        set_usefilter_flag(true);
    }

    void set_overlap_len(unsigned pOverlapLen) { mOverlapLen = pOverlapLen; }

    void set_usefilter_flag(bool pUseFilter) { mUseFilter = pUseFilter; }

    virtual void init_from_args() {}

    // Method prototype
    int filter_sites(mikan::MKSeedSites &pSeedSites, mikan::MKRMAWithSites &pRNAWithSites,
                     mikan::MKSiteScores &pSiteScores);

protected:
    // Define method
    virtual float get_precedence(unsigned, mikan::MKSeedSites &, mikan::MKSiteScores &) { return 0; }

    // Method prototypes
    virtual void set_intervals(mikan::MKSeedSites &pSeedSites, mikan::MKSiteScores &pSiteScores,
                               unsigned pSiteIdx, unsigned &pStartSearch, unsigned &pEndSearch,
                               unsigned &pStartAdd, unsigned &pEndAdd, bool &pSearchOverlap);

    virtual void mark_sites(mikan::MKSeedSites &pSeedSites, mikan::TMRNAPosSet &pSortedPos,
                            mikan::MKSiteScores &pSiteScores);

    void sort_sites(mikan::TMRNAPosSet &pSortedPos, mikan::MKSeedSites &pSeedSites,
                    mikan::MKSiteScores &pSiteScores, mikan::TMRNAPosSet &pSortedSites);

    // Define types
    typedef std::multimap<float, unsigned> TPosMap;
    typedef std::multimap<float, unsigned>::iterator TItMap;
    typedef std::pair<float, unsigned> TPosPair;

    // Define variables
    unsigned mOverlapLen;
    bool mUseFilter;
    mikan::MKOptions const &mOpts;

};

//
// Filter sites by N top scored sites
//
class MKTopNSites : public MKSiteFilter {
public:
    // Define methods
    explicit MKTopNSites(mikan::MKOptions const &opts) : MKSiteFilter(opts), mTopN(0) {}

    virtual void init_from_args();

protected:
    // Define method
    virtual float get_precedence(unsigned pSitePos, mikan::MKSeedSites &pSeedSites,
                                 mikan::MKSiteScores &pSiteScores);

private:
    int mTopN;

    virtual void mark_sites(mikan::MKSeedSites &pSeedSites, mikan::TMRNAPosSet &pSortedPos,
                    mikan::MKSiteScores &pSiteScores);


};

} // namespace mikan

#endif /* MK_SITE_FILTER_HPP_ */
