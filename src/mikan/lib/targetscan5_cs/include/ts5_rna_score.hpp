#ifndef TS5_RNA_SCORE_HPP_
#define TS5_RNA_SCORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "ts5_site_score.hpp"          // TS5SiteScores

namespace ts5cs {

//
// Total context scores
//
class TS5TotalScores {
public:
    // Define methods
    TS5TotalScores() {}

    const seqan::String<float> &get_scores() { return mTotalScores; }

    const seqan::String<int> &get_mrna_pos() { return mMRNAPos; }

    const seqan::String<int> &get_site_num() { return mSiteNum; }

    // Method prototypes
    void clear_scores();

    int calc_scores(TS5SeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKRMAWithSites &pRNAWithSites, TS5SiteScores &pSiteScores);

private:
    seqan::String<float> mTotalScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace ts5cs

#endif /* TS5_RNA_SCORE_HPP_ */
