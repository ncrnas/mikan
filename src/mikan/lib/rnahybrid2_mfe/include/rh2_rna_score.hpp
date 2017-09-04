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
    explicit RH2RNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts) {
        resize(mScoreTypes, 3);
        mScoreTypes[0] = "mfe";
        mScoreTypes[1] = "tfe";
        mScoreTypes[2] = "tnr";
    }

    const seqan::String<float> &get_norm_scores() { return mNormScores; }

    // Method prototypes
    virtual void clear_scores();

    virtual float get_score(int pTypeIdx, int pIdx) {
        if (pTypeIdx == 1) {
            return mRNAScores[pIdx];
        } else if (pTypeIdx == 2) {
            return mNormScores[pIdx];
        }

        return mMinScores[pIdx];
    }

    using mikan::MKRNAScores::calc_scores;

    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                            mikan::MKRMAWithSites &pRNAWithSites, RH2SiteScores &pMFEScores);

private:
    seqan::String<float> mNormScores;
    seqan::String<float> mMinScores;
    seqan::String<float> mMaxNormScores;

};

} // namespace rh2mfe

#endif /* RH2_RNA_SCORE_HPP_ */
