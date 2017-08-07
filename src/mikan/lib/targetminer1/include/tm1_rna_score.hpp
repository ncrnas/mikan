#ifndef TM1_RNA_SCORE_HPP_
#define TM1_RNA_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_score.hpp"       // MKRNAScores
#include "tm1_mrna_svm.hpp"       // TM1MRNAInputVector
#include "tm1_site_filter.hpp"    // TM1SortedSitePos

namespace tm1p {

//
// Classified scores
//
class TM1RNAScores : public mikan::MKRNAScores {
public:
    // Define methods
    explicit TM1RNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts) {}

    const seqan::String<int> &get_labels() { return mPredictions; }

    // Method prototypes
    virtual void clear_scores();

    using mikan::MKRNAScores::calc_scores;
    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKRMAWithSites &pRNAWithSites, TM1SiteScores &mSiteScores);

private:
    seqan::String<int> mPredictions;

    TM1MRNAFeatures mMRNAFeatures;
    TM1MRNAInputVector mMRNAInput;

};

} // namespace tm1p

#endif /* TM1_RNA_SCORE_HPP_ */
