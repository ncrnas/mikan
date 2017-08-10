#ifndef MK_CORE_TMPL_HPP_
#define MK_CORE_TMPL_HPP_

#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "mk_sequence.hpp"        // MKSequences
#include "mk_option.hpp"          // MKOptions
#include "mk_core_base.hpp"       // MKCoreBase
#include "mk_rna_sites.hpp"       // MKRMAWithSites
#include "mk_site_filter.hpp"     // MKTopNSites
#include "mk_rna_score.hpp"       // MKRNAScores

namespace mikan {

//
// MK score process core template
//
template<class TSeedSeqs, class TSeedSites, class TSiteScores, class TSiteFilter, class TRNAScores>
class MKCoreTmpl : public mikan::MKCoreBase {
public:
    // Define methods
    MKCoreTmpl(mikan::MKOptions const &pOpts,
               mikan::TCharSet const &pMiRNAIds,
               mikan::TRNASet const &pMiRNASeqs,
               mikan::TCharSet const &pMRNAIds,
               mikan::TRNASet const &pMRNASeqs,
               mikan::TIndexQGram &pRNAIdx,
               mikan::TFinder &pFinder) :
            MKCoreBase(pOpts, pMiRNAIds, pMiRNASeqs, pMRNAIds, pMRNASeqs, pRNAIdx, pFinder),
            mSeedSeqs(pOpts),
            mSeedSites(pRNAIdx, pFinder, pMRNASeqs),
            mSiteScores(pOpts),
            mRNAWithSites(pOpts),
            mSiteFilter(pOpts),
            mTopNSites(pOpts),
            mRNAScores(pOpts) {}

    virtual int find_seed_sites(unsigned pIdx) {
        int retVal;
        mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

        if (mFindSeedSites) {
            retVal = mSeedSeqs.create_seed_seqs(miRNASeq);
            if (retVal != 0) {
                return 1;
            }

            retVal = mSeedSites.find_seed_sites(mSeedSeqs);
            if (retVal != 0) {
                return 1;
            }

        }

        if (mFilterSites) {
            if (mClusterSites1) {
                mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
            }

            retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
            if (retVal != 0) {
                return 1;
            }
        }

        return 0;
    }

    virtual int calc_site_scores(unsigned pIdx) {
        int retVal;
        mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

        if (mCalcSiteScore) {
            retVal = mSiteScores.calc_scores(miRNASeq, mMRNASeqs, mSeedSites, mRNAWithSites);
            if (retVal != 0) {
                return 1;
            }
        }

        if (mFilterSiteScores) {
            if (mClusterSites2) {
                mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
            }

            retVal = mSiteFilter.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
            if (retVal != 0) {
                return 1;
            }

            if (mSelectTopSites) {
                retVal = mTopNSites.filter_sites(mSeedSites, mRNAWithSites, mSiteScores);
                if (retVal != 0) {
                    return 1;
                }
            }
        }

        return 0;
    }

    virtual int calc_rna_scores(unsigned) {
        int retVal;

        if (mCalcRNAScore) {
            if (mClusterSites3) {
                mRNAWithSites.create_mrna_site_map(mSeedSites, mSiteScores);
            }

            retVal = mRNAScores.calc_scores(mSeedSites, mMRNASeqs, mRNAWithSites, mSiteScores);
            if (retVal != 0) {
                return 1;
            }
        }

        return 0;
    }

    virtual int output_results(unsigned pIdx) {
        int retVal;
        mikan::TRNAStr miRNASeq = mMiRNASeqs[pIdx];

        // Write site scores
        if (mOutputSite) {
            retVal = write_site_score(mMiRNAIds[pIdx]);
            if (retVal != 0) { ;
                return 1;
            }
        }

        // Write total scores
        if (mOutputRNA) {
            retVal = write_rna_score(mMiRNAIds[pIdx]);
            if (retVal != 0) {
                return 1;
            }
        }

        // Write alignments
        if (mOutputAlign) {
            retVal = write_alignment(mMiRNAIds[pIdx]);
            if (retVal != 0) {
                return 1;
            }
        }

        return 0;
    }

    virtual void clear_all() {
        mSeedSeqs.clear_seeds();
        mSeedSites.clear_pos();
        mSiteScores.clear_scores();
        mRNAWithSites.clear_maps();
        mRNAScores.clear_scores();
    }

protected:
    TSeedSeqs mSeedSeqs;
    TSeedSites mSeedSites;
    TSiteScores mSiteScores;
    mikan::MKRMAWithSites mRNAWithSites;
    TSiteFilter mSiteFilter;
    mikan::MKTopNSites mTopNSites;
    TRNAScores mRNAScores;

    virtual int write_site_score(seqan::CharString const &pMiRNAId) = 0;

    virtual int write_rna_score(seqan::CharString const &pMiRNAId) = 0;

    virtual int write_alignment(seqan::CharString const &pMiRNAId) = 0;

};

} // namespace mikan

#endif /* MK_CORE_TMPL_HPP_ */
