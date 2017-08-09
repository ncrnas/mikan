#ifndef TS5_CORE_HPP_
#define TS5_CORE_HPP_

#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_core.hpp"            // MKCoreBase
#include "ts5_feature.hpp"        // TS5RawFeatures
#include "ts5_option.hpp"         // TS5CSOptions
#include "ts5_site_score.hpp"     // TS5SiteScores
#include "ts5_seed_site.hpp"      // TS5SeedSites
#include "ts5_rna_score.hpp"      // TS5RNAScores

namespace ts5cs {

int TS5CoreMain(int argc, char const **argv);

//
// TargetScan context score process core
//
class TS5Core : public mikan::MKCoreBase {
public:
    // Define methods
    TS5Core(mikan::MKOptions const &pOpts, mikan::TCharSet const &pMiRNAIds, mikan::TRNASet const &pMiRNASeqs,
    mikan::TCharSet const &pMRNAIds, mikan::TRNASet const &pMRNASeqs,
    mikan::TIndexQGram &pRNAIdx, mikan::TFinder &pFinder) :
    MKCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
    mSeedSeqs(pOpts), mSeedSites(pRNAIdx, pFinder, pMRNASeqs), mSiteScores(pOpts),
    mRNAWithSites(pOpts), mRNAScores(pOpts) {

        mFilterSites = false;
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
    TS5SeedSeqs mSeedSeqs;
    TS5SeedSites mSeedSites;
    TS5SiteScores mSiteScores;
    mikan::MKRMAWithSites mRNAWithSites;
    TS5RNAScores mRNAScores;

private:
    int write_site_score(seqan::CharString const &pMiRNAId);
    int write_rna_score(seqan::CharString const &pMiRNAId);
    int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace ts5cs

#endif /* TS5_CORE_HPP_ */
