#ifndef MKE_RNA_SCORE_HPP_
#define MKE_RNA_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_score.hpp"       // MKRNAScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_rna_score.hpp"       // MKRNAScores
#include "mke_site_score.hpp"     // MKESiteScores
#include "mke_rna_sites.hpp"      // MKERMAWithSites

namespace mkens {

//
// Total context scores
//
class MKERNAScores : public mikan::MKRNAScores {
public:
    // Define methods
    explicit MKERNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts) {
        mScoreTypeN = 0;
    }

    // Method prototypes
    virtual void clear_scores();

    seqan::CharString &get_tool_score(int pIdx) { return mToolScores[pIdx]; }

    void add_score_types(mkens::MKEOptions const &pMKEOpts, mikan::MKRNAScores &pRNAScores,
                         seqan::CharString &pPrefix);

    void init_score_list(mkens::MKERMAWithSites &pRNAWithSites);

    void add_scores(MKEOptions const &pMKEOpts, mkens::MKERMAWithSites &pRNAWithSites,
                    mikan::MKRNAScores &pRNAScores, seqan::CharString &pPrefix);

    void combine_scores(MKEOptions const &pMKEOpts);

    void set_site_count(mikan::MKSeedSites &pSeedSites, mikan::MKSiteScores &pSiteScores,
                        mkens::MKERMAWithSites &pRNAWithSites);

    // Variables
    unsigned mScoreTypeN;
    std::map<std::string, unsigned> mIdxMap;

    seqan::StringSet<seqan::StringSet<float> > mRNARawScoreList;
    seqan::StringSet<seqan::StringSet<float> > mRNANormScoreList;
    mikan::TCharSet mToolScores;

private:
    float normalize_score(float pScore, MKEOptions const &pMKEOpts, seqan::CharString &pScoreType);

};

} // namespace mkens

#endif /* MKE_RNA_SCORE_HPP_ */
