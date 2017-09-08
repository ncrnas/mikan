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
        resize(mScoreTypes, 2);
        mScoreTypes[0] = "mfe";
        mScoreTypes[1] = "nrm";
    }

    const seqan::String<float> &get_norm_scores() { return mNormScores; }

    const seqan::String<float> &get_mfe_minlogtotal() { return mLogMinRNAScores; }

    const seqan::String<float> &get_norm_maxlogtotal() { return mLogMaxNormScores; }

    // Method prototypes
    virtual void clear_scores();

    virtual float get_score(int pTypeIdx, int pIdx) {
        if (pTypeIdx == 1) {
            return mLogMaxNormScores[pIdx];
        }

        return mLogMinRNAScores[pIdx];
    }

    using mikan::MKRNAScores::calc_scores;

    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                            mikan::MKRMAWithSites &pRNAWithSites, RH2SiteScores &pMFEScores);

private:
    seqan::String<float> mNormScores;
    seqan::String<float> mMinScores;
    seqan::String<float> mMaxNormScores;
    seqan::String<float> mLogMinRNAScores;
    seqan::String<float> mLogMaxNormScores;

};

} // namespace rh2mfe

#endif /* RH2_RNA_SCORE_HPP_ */
