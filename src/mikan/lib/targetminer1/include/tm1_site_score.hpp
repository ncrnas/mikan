#ifndef TM1_SITE_SCORE_HPP_
#define TM1_SITE_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_option.hpp"          // MKOptions
#include "mk_seed_site.hpp"       // MKSeedSites
#include "tm1_site_feature.hpp"   // TM1SiteFeatures

namespace tm1p {

//
// Store site level scores
//
class TM1SiteScores : public mikan::MKSiteScores {
public:
    // Define methods
    TM1SiteScores(mikan::MKOptions const &opts) : MKSiteScores(opts) {}

    void write_alignment(int pIdx) { mSiteFeatures.write_alignment(pIdx); }

    // Method prototypes
    void clear_scores();

    int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites, mikan::MKRMAWithSites &pRNAWithSites);

    TM1SiteFeatures &get_site_features() { return mSiteFeatures; }

private:
    TM1SiteFeatures mSiteFeatures;

};

} // namespace tm1p

#endif /* TM1_SITE_SCORE_HPP_ */
