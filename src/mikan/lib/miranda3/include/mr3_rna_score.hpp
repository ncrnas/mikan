#ifndef MR3_RNA_SCORE_HPP_
#define MR3_RNA_SCORE_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mr3_option.hpp"         // MR3Options
#include "mr3_site_score.hpp"     // MR3SiteScores
#include "mr3_seed_site.hpp"      // MR3SeedSites

namespace mr3as {

//
// Total scores
//
class MR3TotalScores {
public:
    // Constant values
    const double MIN_EXP_DIFF;

public:
    // Define methods
    MR3TotalScores() : MIN_EXP_DIFF(-100.0) {}

    const seqan::String<float> &get_align_scores() { return mTotalAlignScores; }

    const seqan::String<float> &get_energy_scores() { return mTotalEnScores; }

    const seqan::String<int> &get_mrna_pos() { return mMRNAPos; }

    const seqan::String<int> &get_site_num() { return mSiteNum; }

    // Method prototype
    void clear_scores();

    int calc_scores(MR3SeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKRMAWithSites &pRNAWithSites, MR3SiteScores &pSiteScores);

private:
    seqan::String<float> mTotalAlignScores;
    seqan::String<float> mTotalEnScores;
    seqan::String<float> mLogMaxAlignScores;
    seqan::String<float> mLogMaxEnScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace mr3as

#endif /* MR3_RNA_SCORE_HPP_ */
