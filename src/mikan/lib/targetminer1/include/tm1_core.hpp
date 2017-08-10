#ifndef TM1_CORE_HPP_
#define TM1_CORE_HPP_

#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "mk_site_score.hpp"     // MKSiteScores
#include "mk_option.hpp"         // MKOptions
#include "mk_core_tmpl.hpp"      // MKCoreTmpl
#include "tm1_option.hpp"        // TM1Options
#include "tm1_site_filter.hpp"   // TM1SiteFilter
#include "tm1_mrna_feature.hpp"  // TM1MRNAFeatures
#include "tm1_mrna_svm.hpp"      // TM1MRNAModel, TM1MRNAInputVector
#include "tm1_option.hpp"        // TM1Options
#include "tm1_rna_score.hpp"     // TM1RNAScores
#include "tm1_seed_site.hpp"     // TM1SeedSites
#include "tm1_site_score.hpp"    // TM1SiteScores

namespace tm1p {

//
// TargetScan context score process core
//
typedef mikan::MKCoreTmpl<TM1SeedSeqs, TM1SeedSites, TM1SiteScores, TM1SiteFilter, TM1RNAScores> TM1CoreBase;

class TM1Core : public TM1CoreBase {
public:
    // Define methods
    TM1Core(mikan::MKOptions const &pOpts,
            mikan::TCharSet const &pMiRNAIds,
            mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds,
            mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx,
            mikan::TFinder &pFinder) :
            TM1CoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mClusterSites2 = false;
        mFilterSiteScores = false;
        mSelectTopSites = false;
        mClusterSites3 = false;

    }

private:
    virtual int write_site_score(seqan::CharString const &pMiRNAId);

    virtual int write_rna_score(seqan::CharString const &pMiRNAId);

    virtual int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tm1p

#endif /* TM1_CORE_HPP_ */
