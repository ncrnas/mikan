#ifndef MR3_CORE_HPP_
#define MR3_CORE_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_core_tmpl.hpp"       // MKCoreTmpl
#include "mr3_option.hpp"         // MR3Options
#include "mr3_site_score.hpp"     // MR3GGDScores
#include "mr3_seed_site.hpp"      // MR3SeedSites
#include "mr3_site_filter.hpp"    // MR3SiteFilter
#include "mr3_rna_score.hpp"      // MR3RNAScores

namespace mr3as {

int MR3CoreMain(int argc, char const **argv);

//
// MR3 score process core
//
typedef mikan::MKCoreTmpl<MR3SeedSeqs, MR3SeedSites, MR3SiteScores, MR3SiteFilter, MR3RNAScores > MR3CoreBase;

class MR3Core : public MR3CoreBase {
public:
    // Define methods
    MR3Core(mikan::MKOptions const &pOpts,
            mikan::TCharSet const &pMiRNAIds,
            mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds,
            mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx,
            mikan::TFinder &pFinder) :
            MR3CoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mClusterSites1 = false;
        mFilterSites = false;
        mSelectTopSites = false;
        mClusterSites3 = false;

    }

private:
    virtual int write_site_score(seqan::CharString const &pMiRNAId);
    virtual int write_rna_score(seqan::CharString const &pMiRNAId);
    virtual int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace mr3as

#endif /* MR3_CORE_HPP_ */
