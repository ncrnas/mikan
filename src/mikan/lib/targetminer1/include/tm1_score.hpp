#ifndef TM1_SCORE_HPP_
#define TM1_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_typedef.hpp"         // TRNATYPE
#include "tm1_mrna_svm.hpp"       // TM1MRNAInputVector
#include "tm1_site_cluster.hpp"   // TM1SortedSitePos

namespace tm1p {

//
// Classified scores
//
class TM1ClassifiedScores {
public:
    // Define methods
    TM1ClassifiedScores() {}

    const seqan::String<float> &get_scores() { return mScores; }

    const seqan::String<int> &get_labels() { return mPredictions; }

    const seqan::String<unsigned> &get_site_num() { return mSiteNum; }

    // Method prototypes
    void clear_scores();

    int calc_scores(const seqan::String<unsigned> &pSiteCoutns, const seqan::String<float> &pScores);

private:
    seqan::String<float> mScores;
    seqan::String<int> mPredictions;
    seqan::String<unsigned> mSiteNum;
};

} // namespace tm1p

#endif /* TM1_SCORE_HPP_ */
