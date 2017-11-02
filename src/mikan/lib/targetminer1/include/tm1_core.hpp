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

        mFindSeedSites = true;
        mFilterSites = true;
        mCalcSiteScore = true;
        mClusterSites1 = true;
        mFilterSiteScores = false;
        mClusterSites2 = false;
        mSelectTopSites = false;
        mClusterSites3 = false;
        mCalcRNAScore = true;

    }

private:
    virtual void prepare_site_output(mikan::TCharStr const &pMiRNAId, unsigned pRNAPosIdx, unsigned pSitePosIdx);

    virtual void prepare_rna_output(mikan::TCharStr const &pMiRNAId);

    virtual int write_alignment(mikan::TCharStr const &pMiRNAId);

};

} // namespace tm1p

#endif /* TM1_CORE_HPP_ */
