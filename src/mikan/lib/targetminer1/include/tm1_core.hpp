#ifndef TM1_CORE_HPP_
#define TM1_CORE_HPP_

#include "mk_typedef.hpp"        // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"       // MKSequences
#include "mk_site_score.hpp"     // MKSiteScores
#include "mk_option.hpp"         // MKOptions
#include "mk_core.hpp"           // MKCoreBase
#include "tm1_site_filter.hpp"   // TM1SiteFilter
#include "tm1_mrna_feature.hpp"  // TM1MRNAFeatures
#include "tm1_mrna_svm.hpp"      // TM1MRNAModel, TM1MRNAInputVector
#include "tm1_option.hpp"        // TM1CSOptions
#include "tm1_rna_score.hpp"     // TM1RNAScores
#include "tm1_seed_site.hpp"     // TM1SeedSites
#include "tm1_site_score.hpp"    // TM1SiteScores

namespace tm1p {

int TM1CoreMain(int argc, char const **argv);

//
// TargetScan context score process core
//
class TM1Core : public mikan::MKCoreBase {
public:
    // Define methods
    TM1Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
    mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
    mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
    MKCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
    mSeedSeqs(pOpts), mSeedSites(pRNAIdx, pFinder, pMRNASeqs), mRNAWithSites(pOpts),
    mSiteScores(pOpts), mSiteFilter(pOpts), mTopNSites(pOpts), mRNAScores(pOpts) {

        mFilterSiteScores = false;
        mSelectTopSites = false;
        mSelectTopRNAs = false;

    }

    virtual int find_seed_sites(unsigned pIdx);
    virtual int calc_site_scores(unsigned pIdx);
    virtual int ensemble_site_scores(unsigned pIdx);
    virtual int calc_rna_scores(unsigned pIdx);
    virtual int ensemble_rna_scores(unsigned pIdx);
    virtual int output_results(unsigned pIdx);
    virtual void clear_all();

private:
    TM1SeedSeqs mSeedSeqs;
    TM1SeedSites mSeedSites;
    mikan::MKRMAWithSites mRNAWithSites;
    TM1SiteScores mSiteScores;
    TM1SiteFilter mSiteFilter;
    mikan::MKTopNSites mTopNSites;
    TM1RNAScores mRNAScores;

    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_rna_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace tm1p

#endif /* TM1_CORE_HPP_ */
