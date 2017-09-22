#ifndef MK_RNA_SCORE_HPP_
#define MK_RNA_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_option.hpp"          // MKOptions
#include "mk_seed_site.hpp"       // MKSeedSites
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites

namespace mikan {

//
// Store RNA level scores
//

class MKRNAScores {
public:

    // Define variable
    seqan::String<bool> mEffectiveRNAs;

    // Define methods
    explicit MKRNAScores(mikan::MKOptions const &opts) : mOpts(opts) {
        resize(mScoreTypes, 0);
    }

    virtual const seqan::String<float> &get_scores() { return mRNAScores; }

    virtual float get_score(int, int pIdx) { return mRNAScores[pIdx]; }

    const mikan::TMRNAPosSet &get_mrna_pos() { return mMRNAPos; }

    const seqan::String<unsigned> &get_site_num() { return mSiteNum; }

    // Method prototypes
    virtual void clear_scores();

    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                            mikan::MKRMAWithSites &pRNAWithSites, mikan::MKSiteScores &pSiteScores);

    const mikan::TCharSet &get_score_types() { return mScoreTypes; }

protected:
    // Define variables
    mikan::MKOptions const &mOpts;
    seqan::String<float> mRNAScores;
    mikan::TMRNAPosSet mMRNAPos;
    seqan::String<unsigned> mSiteNum;
    mikan::TCharSet mScoreTypes;

};

} // namespace mikan

#endif /* MK_RNA_SCORE_HPP_ */
