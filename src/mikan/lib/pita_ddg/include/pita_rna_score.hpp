#ifndef PITA_RNA_SCORE_HPP_
#define PITA_RNA_SCORE_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_seed_site.hpp"       // MKSeedSites
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_rna_score.hpp"       // MKRNAScores
#include "pita_option.hpp"        // PITAOptions
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_score.hpp"    // PITAGGDScores
#include "vr16_ddg_core.hpp"      // VR16DDGWorkSpace

namespace ptddg {

//
// Total ddG scores
//
class PITARNAScores : public mikan::MKRNAScores {
public:
    // Constant values
    const double MIN_EXP_DIFF;

    // Define type
    typedef std::set<unsigned>::iterator TItSet;

    // Define methods
    explicit PITARNAScores(mikan::MKOptions const &opts) : MKRNAScores(opts), MIN_EXP_DIFF(-100.0) {}

    // Method prototype
    using mikan::MKRNAScores::calc_scores;

    virtual int calc_scores(mikan::MKSeedSites &pSeedSites, mikan::TRNASet const &pMRNASeqs,
                            mikan::MKRMAWithSites &pRNAWithSites, PITASiteScores &pDDGScores);

};

} // namespace ptddg

#endif /* PITA_RNA_SCORE_HPP_ */
