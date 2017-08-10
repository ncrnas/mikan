#ifndef RH2_RNA_SCORE_HPP_
#define RH2_RNA_SCORE_HPP_

#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_rna_score.hpp"       // MKRNAScores
#include "rh2_site_score.hpp"     // RH2SiteScores

namespace rh2mfe {

//
// Total MFE scores
//
class RH2RNAScores : public mikan::MKRNAScores {
public:
    // Define methods
    explicit RH2RNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts) {}

    const seqan::String<float> &get_norm_scores() { return mNormScores; }

    // Method prototypes
    virtual void clear_scores();

    using mikan::MKRNAScores::calc_scores;

    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                            mikan::MKRMAWithSites &pRNAWithSites, RH2SiteScores &pMFEScores);

private:
    seqan::String<float> mNormScores;

};

} // namespace rh2mfe

#endif /* RH2_RNA_SCORE_HPP_ */
