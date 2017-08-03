#ifndef PITA_RNA_SCORE_HPP_
#define PITA_RNA_SCORE_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "pita_option.hpp"        // PITAOptions
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_score.hpp"    // PITAGGDScores
#include "vr16_ddg_core.hpp"      // VR16DDGWorkSpace

namespace ptddg {

//
// Total ddG scores
//
class PITATotalScores {
public:
    // Constant values
    const double MIN_EXP_DIFF;

    typedef std::set<unsigned>::iterator TItSet;

public:
    // Define methods
    PITATotalScores() : MIN_EXP_DIFF(-100.0) {}

    const seqan::String<float> &get_scores() { return mTotalScores; }

    const seqan::String<int> &get_mrna_pos() { return mMRNAPos; }

    const seqan::String<int> &get_site_num() { return mSiteNum; }

    // Method prototype
    void clear_scores();

    int calc_scores(PITASeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                    mikan::MKRMAWithSites &pRNAWithSites, PITASiteScores &pDDGScores);

private:
    seqan::String<float> mTotalScores;
    seqan::String<int> mMRNAPos;
    seqan::String<int> mSiteNum;

};

} // namespace ptddg

#endif /* PITA_RNA_SCORE_HPP_ */
