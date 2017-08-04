#ifndef MK_SITE_SCORE_HPP_
#define MK_SITE_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_option.hpp"          // MKOptions
#include "mk_seed_site.hpp"       // MKSeedSites

namespace mikan {

//
// Store site level scores
//
class MKRMAWithSites;  // To avoid a cross-reference issue between MKSiteScores and MKRMAWithSites

class MKSiteScores {
public:
    // Constant value
    static const unsigned INDEXED_SEQ_LEN = mikan::SEEDLEN;

    // Define variable
    seqan::String<bool> mEffectiveSites;

    // Define methods
    explicit MKSiteScores(mikan::MKOptions const &opts) : mOpts(opts) {}

    virtual float get_score(int pIdx) { return mSiteScores[pIdx]; }

    virtual void init_from_args() {}

    virtual int get_wide_site_start(int) { return 0; }

    virtual int get_wide_site_length(int) { return 0; }

    // Method prototypes
    virtual void clear_scores();

    virtual int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites, mikan::MKRMAWithSites const &pRNAWithSites);

protected:
    // Define variables
    mikan::MKOptions const &mOpts;
    mikan::TScoreSet mSiteScores;

};

} // namespace mikan

#endif /* MK_SITE_SCORE_HPP_ */
