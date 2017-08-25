#ifndef MKE_SITE_SCORE_HPP_
#define MKE_SITE_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mke_seed_site.hpp"      // MKESeedSites

namespace mkens {

//
// Store site scores
//
class MKESiteScores : public mikan::MKSiteScores {
public:
    // Define methods
    explicit MKESiteScores(mikan::MKOptions const &opts) :
            MKSiteScores(opts) {}

};

} // namespace mkens

#endif /* MKE_SITE_SCORE_HPP_ */

