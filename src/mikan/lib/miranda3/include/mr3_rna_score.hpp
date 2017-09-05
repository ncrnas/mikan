#ifndef MR3_RNA_SCORE_HPP_
#define MR3_RNA_SCORE_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_site.hpp"       // MKSeedSites
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_rna_score.hpp"       // MKRNAScores
#include "mr3_option.hpp"         // MR3Options
#include "mr3_site_score.hpp"     // MR3SiteScores

namespace mr3as {

//
// Total scores
//
class MR3RNAScores : public mikan::MKRNAScores {
public:
    // Constant values
    const double MIN_EXP_DIFF;

public:
    // Define methods
    explicit MR3RNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts), MIN_EXP_DIFF(-100.0) {
        resize(mScoreTypes, 2);
        mScoreTypes[0] = "alg";
        mScoreTypes[1] = "eng";
    }

    const seqan::String<float> &get_align_scores() { return mTotalAlignScores; }

    const seqan::String<float> &get_energy_scores() { return mTotalEnScores; }

    const seqan::String<float> &get_align_maxlogtotal() { return mLogMaxAlignScores; }

    const seqan::String<float> &get_energy_minlogtotal() { return mLogMinEnScores; }

    virtual float get_score(int pTypeIdx, int pIdx) {
        if (pTypeIdx == 1) {
            return mLogMaxAlignScores[pIdx];
        }

        return mLogMinEnScores[pIdx];
    }

    // Method prototype
    virtual void clear_scores();

    using mikan::MKRNAScores::calc_scores;

    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                            mikan::MKRMAWithSites &pRNAWithSites, MR3SiteScores &pSiteScores);

private:
    seqan::String<float> mTotalAlignScores;
    seqan::String<float> mTotalEnScores;
    seqan::String<float> mLogMaxAlignScores;
    seqan::String<float> mLogMinEnScores;
    seqan::String<float> mMaxAlignScores;
    seqan::String<float> mMinEnScores;

};

} // namespace mr3as

#endif /* MR3_RNA_SCORE_HPP_ */
