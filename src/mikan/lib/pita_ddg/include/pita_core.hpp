#ifndef PITA_CORE_HPP_
#define PITA_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_core_tmpl.hpp"       // MKCoreTmpl
#include "pita_option.hpp"        // PITAOptions
#include "pita_site_score.hpp"    // PITAGGDScores
#include "pita_seed_site.hpp"     // PITASeedSites
#include "pita_site_filter.hpp"   // PITASiteFilter
#include "pita_rna_score.hpp"     // PITARNAScores

namespace ptddg {

//
// PITA score process core
//
typedef mikan::MKCoreTmpl<PITASeedSeqs, PITASeedSites, PITASiteScores, PITASiteFilter, PITARNAScores> PITACoreBase;

class PITACore : public PITACoreBase {
public:
    // Define methods
    PITACore(mikan::MKOptions const &pOpts,
             mikan::TCharSet const &pMiRNAIds,
             mikan::TRNASet const &pMiRNASeqs,
             mikan::TCharSet const &pMRNAIds,
             mikan::TRNASet const &pMRNASeqs,
             mikan::TIndexQGram &pRNAIdx,
             mikan::TFinder &pFinder) :
            PITACoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mClusterSites2 = false;
        mFilterSiteScores = false;
        mSelectTopSites = false;
        mClusterSites3 = false;

    }

private:
    virtual int write_site_score(mikan::TCharStr const &pMiRNAId);

    virtual int write_site_score_gff(mikan::TCharStr const &) { return 0; }

    virtual int write_rna_score(mikan::TCharStr const &pMiRNAId);

    virtual int write_rna_score_gff(mikan::TCharStr const &) { return 0; }

    virtual int write_alignment(mikan::TCharStr const &pMiRNAId);

};

} // namespace ptddg

#endif /* PITA_CORE_HPP_ */
