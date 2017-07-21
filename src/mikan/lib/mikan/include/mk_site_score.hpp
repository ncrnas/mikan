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
class MKSiteScores {
public:
    // Constant value
    static const unsigned INDEXED_SEQ_LEN = mikan::SEEDLEN;

    // Define variable
    seqan::String<bool> mEffectiveSites;

    // Define methods
    explicit MKSiteScores(mikan::MKOptions const &opts) :
            mOpts(opts) {}

    virtual float get_score(int pIdx) { return mSiteScores[pIdx]; }

    void init_from_args() {}

    // Method prototypes
    void clear_scores();

    int calc_scores(mikan::TRNAStr const &pMiRNASeq, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKSeedSites &pSeedSites);

protected:
    // Define variable
    mikan::MKOptions const &mOpts;

    mikan::TScoreSet mSiteScores;

};

} // namespace mikan

#endif /* MK_SITE_SCORE_HPP_ */
