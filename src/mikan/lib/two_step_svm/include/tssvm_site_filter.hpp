#ifndef TSSVM_SITE_FILTER_HPP_
#define TSSVM_SITE_FILTER_HPP_

#include <seqan/sequence.h>
#include <seqan/index.h>
#include "mk_typedef.hpp"           // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_seq.hpp"          // MKSeedSeqs
#include "mk_seed_site.hpp"         // MKSeedSites
#include "mk_site_filter.hpp"       // MKSiteFilter
#include "mk_option.hpp"            // MKOptions

namespace tssvm {

//
// miRNA seed sites overlap
//
class TSSVMSiteFilter : public mikan::MKSiteFilter {
public:
    // Define types
    typedef std::set<unsigned>::iterator TItSet;
    typedef std::multimap<unsigned, unsigned>::iterator TItMap;
    typedef std::pair<unsigned, unsigned> TPosPair;

    // Define methods
    TSSVMSiteFilter(mikan::MKOptions const &opts) : MKSiteFilter(opts) {}

private:
    float get_precedence(unsigned pSitePos, mikan::MKSeedSites &pSeedSites,
                         mikan::MKSiteScores &pSiteScores);

    void set_intervals(mikan::MKSeedSites &pSeedSites, mikan::MKSiteScores &pSiteScores,  unsigned pSiteIdx,
                       unsigned &pStartSearch, unsigned &pEndSearch, unsigned &pStartAdd, unsigned &pEndAdd,
                       bool &pSearchOverlap);

};

} // namespace tssvm

#endif /* TSSVM_SITE_FILTER_HPP_ */
