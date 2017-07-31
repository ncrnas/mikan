#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mr3_score.hpp"          // MR3MFEScores
#include "mr3_site_filter.hpp"    // MR3SiteFilter

using namespace seqan;

namespace mr3as {

//
// MR3SiteFilter methods
//
float MR3SiteFilter::get_precedence(
        unsigned pSitePos,
        mikan::MKSeedSites &,
        mikan::MKSiteScores &pSiteScores) {

    float preced = pSiteScores.get_score(pSitePos) * -1.0;

    return preced;
}

} // namespace mr3as
