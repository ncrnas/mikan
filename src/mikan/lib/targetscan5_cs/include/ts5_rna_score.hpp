#ifndef TS5_RNA_SCORE_HPP_
#define TS5_RNA_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_score.hpp"       // MKRNAScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_rna_score.hpp"       // MKRNAScores
#include "ts5_site_score.hpp"     // TS5SiteScores

namespace ts5cs {

//
// Total context scores
//
class TS5RNAScores : public mikan::MKRNAScores {
public:
    // Define methods
    TS5RNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts) {}

    // Method prototypes
    using mikan::MKRNAScores::calc_scores;
    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKRMAWithSites &pRNAWithSites, TS5SiteScores &pDDGScores);

};

} // namespace ts5cs

#endif /* TS5_RNA_SCORE_HPP_ */
