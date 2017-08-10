#ifndef MKE_CORE_HPP_
#define MKE_CORE_HPP_

#include "mk_typedef.hpp"         // TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_site_filter.hpp"     // MKSiteFilter
#include "mk_site_score.hpp"      // MKSiteScores
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_core_tmpl.hpp"       // MKCoreTmpl
#include "mke_option.hpp"         // MKEOptions
#include "mke_site_score.hpp"     // MKESiteScores
#include "mke_seed_site.hpp"      // MKESeedSites
#include "mke_rna_score.hpp"      // MKERNAScores

namespace mkens {

//
// mikan ensemble score process core
//
typedef mikan::MKCoreTmpl<MKESeedSeqs, MKESeedSites, MKESiteScores, mikan::MKSiteFilter, MKERNAScores> MKECoreBase;

class MKECore : public MKECoreBase {
public:
    // Define methods
    MKECore(mikan::MKOptions const &pOpts,
            mikan::TCharSet const &pMiRNAIds,
            mikan::TRNASet const &pMiRNASeqs,
            mikan::TCharSet const &pMRNAIds,
            mikan::TRNASet const &pMRNASeqs,
            mikan::TIndexQGram &pRNAIdx,
            mikan::TFinder &pFinder) :
            MKECoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder) {

        mClusterSites1 = false;
        mFilterSites = false;
        mClusterSites2 = false;
        mFilterSiteScores = false;
        mSelectTopSites = false;

    }

private:
    virtual int write_site_score(seqan::CharString const &pMiRNAId);

    virtual int write_rna_score(seqan::CharString const &pMiRNAId);

    virtual int write_alignment(seqan::CharString const &pMiRNAId);

};

} // namespace mkens

#endif /* MKE_CORE_HPP_ */
