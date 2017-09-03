#ifndef MKE_SITE_SCORE_HPP_
#define MKE_SITE_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mke_seed_site.hpp"      // MKESeedSites

namespace mkens {

//
// Store site scores
//
class MKESiteScores : public mikan::MKSiteScores {
public:
    // Define method
    explicit MKESiteScores(mikan::MKOptions const &pOpts) : MKSiteScores(pOpts) {
        mScoreTypeN = 0;
    }

    // Method prototypes
    virtual void clear_scores();

    seqan::CharString &get_tool_score(int pIdx) { return mToolScores[pIdx]; }

    void add_score_types(mkens::MKEOptions const &pMKEOpts, mikan::MKSiteScores &pSiteScores,
                         seqan::CharString &pPrefix);

    void init_score_list(MKESeedSites &pMKESeedSites);

    void add_scores(MKEOptions const &pMKEOpts, mikan::MKSeedSites &pSeedSites, MKESeedSites &pMKESeedSites,
                    mikan::MKSiteScores &pSeedScores, seqan::CharString &pPrefix);

    void combine_scores(MKEOptions const &pMKEOpts);

private:
    // Variables
    unsigned mScoreTypeN;
    std::map<std::string, unsigned> mIdxMap;

    seqan::StringSet<seqan::StringSet<float> > mSiteRawScoreList;
    seqan::StringSet<seqan::StringSet<float> > mSiteNormScoreList;
    mikan::TCharSet mToolScores;

    float normalize_score(float pScore, MKEOptions const &pMKEOpts, seqan::CharString &pScoreType);

};

} // namespace mkens

#endif /* MKE_SITE_SCORE_HPP_ */

